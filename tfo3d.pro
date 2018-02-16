Include "tfo3d_data.pro";

Dir="res/";
ExtGmsh = ".pos";
ExtGnuplot  = ".dat";
po = "Output 3D/";

Group{

  Core = #CORE;
  Air = #AIR;

  Primary_Helix    = #PRIMARY;
  Secondary0_Helix = #SECONDARY0;
  Secondary1_Helix = #SECONDARY1;

  Primary_In    = #{(PRIMARY+1)};
  Secondary0_In = #{(SECONDARY0+1)};
  Secondary1_In = #{(SECONDARY1+1)};

  Primary_Out    = #{(PRIMARY+2)};
  Secondary0_Out = #{(SECONDARY0+2)};
  Secondary1_Out = #{(SECONDARY1+2)};

  Surf_P_In  = #{IN_PRI};
  Surf_S0_In = #{IN_SEC0};
  Surf_S1_In = #{IN_SEC1};

  Surf_P_Out  = #{OUT_PRI};
  Surf_S0_Out = #{OUT_SEC0};
  Surf_S1_Out = #{OUT_SEC1};


  Primary     = Region[{Primary_Helix, Primary_In, Primary_Out}];
  Secondary0  = Region[{Secondary0_Helix, Secondary0_In, Secondary0_Out}];
  Secondary1  = Region[{Secondary1_Helix, Secondary1_In, Secondary1_Out}];
  Secondary   = Region[{Secondary0,Secondary1}];

  Surf_In  = Region[{Surf_P_In,  Surf_S0_In,  Surf_S1_In}];
  Surf_Out = Region[{Surf_P_Out, Surf_S0_Out, Surf_S1_Out}];

  Winding     = Region[{Primary, Secondary0, Secondary1}];

  DomainC = Region[{}]; // No conduction domain
  Surf_Elec = Region[{}]; // This could be a surface, in or out
  SkinDomainC = Region[{}];

  DomainCC = Region[ {Air, Core, Winding} ];
  DomainS = Region[ {Winding} ];


  Surf_bn0 = #SURF_AIROUT;
  Surf_FixedMVP = Region[{Surf_bn0}];

  Domain = Region[ {DomainC,DomainCC} ];
  DomainDummy = Region[ 123456789 ];
}

Function{
  // RMS current
  IA_pri = 10;
  IA_sec = 0;

  DefineConstant[
    Freq = { 50., Min 0, Max 1e3, Step 1,
      Name "Input/21Frequency [Hz]", Highlight "AliceBlue"},
    Irms_pri = { IA_pri, Min 1, Max 4*IA_pri, Step 2,
      Name "Input/4Coil Parameters/0Primary Current (rms) [A]", Highlight "AliceBlue"},
    Irms_sec = { IA_sec, Min 1, Max 4*IA_sec, Step 2,
      Name "Input/4Coil Parameters/1Secondary Current (rms) [A]", Highlight "AliceBlue"}
    // NbWires = { Nw,
    // Name "Input/4Coil Parameters/1Number of turns", Highlight "AliceBlue"}
    SymmetryFactor = 1
  ];

  radius[Primary]   = rp;
  radius[Secondary] = rs;
  Irms[Primary]   =  Irms_pri;
  Irms[Secondary] = -Irms_sec; //current in different direction
  Ipk[] = Irms[]*Sqrt[2];

  Ap = -(interwire_pri+rp*2+thick_insul*2)*Np/(10*Np-1);
  As =  (interwire_sec+rs*2+thick_insul*2)*Ns/(12*Ns-1);

  //vDir[Region[{Air,Core}]] = Vector [ 0, 0, 0];
  vDir[Primary_Helix] =    Unit [ Vector [ Z[], -Ap ,-X[] ] ] ;
  vDir[Secondary0_Helix] = Unit [ Vector [ Z[], -As ,-X[] ] ] ;
  vDir[Secondary1_Helix] = Unit [ Vector [ Z[], -As ,-X[] ] ] ;

  vDir[Region[{Primary_In,  Secondary0_In,  Secondary1_In}]]  = Vector [0., 0.,-1.] ;
  vDir[Region[{Primary_Out, Secondary0_Out, Secondary1_Out}]] = Vector [1., 0., 0.] ;

  // SurfCoil[] = Pi*radius[]*radius[]; // Ok but better to take the actual section of one of the in or out surfaces
  SurfCoil[Primary]    = SurfaceArea[]{IN_PRI};
  SurfCoil[Secondary0] = SurfaceArea[]{IN_SEC0};
  SurfCoil[Secondary1] = SurfaceArea[]{IN_SEC1};


  js0[] = Ipk[]/SurfCoil[]*vDir[];

  // Material properties

  mu0 = 4.e-7 * Pi ;

  // Use if linear problem
  DefineConstant[
    sigma_coil = { sigma_cu, Label "Conductivity [S/m]",
      Name "Input/4Coil Parameters/5Conductivity", Highlight "AliceBlue"},
    mur_fe = { 2000, Min 100, Max 2000, Step 100,
      Name "Input/42Core relative permeability", Highlight "AliceBlue"}
  ]; //N27 initial permeability = 2000

  nu [ Region[{Air, Winding}] ]  = 1./mu0 ;
  nu [Core]  = 1/(mur_fe*mu0) ;

  sigma[Winding] = sigma_coil ;
  rho[] = 1/sigma[] ;
}


Jacobian {
  { Name Vol ;
    Case {
      { Region All ;       Jacobian Vol ; }
    }
  }
  { Name Sur ;
    Case {
      { Region All ; Jacobian Sur ; }
    }
  }
}

Integration {
  { Name II ;
    Case {
      {
	Type Gauss ;
	Case {
	  { GeoElement Triangle    ; NumberOfPoints  4 ; }
	  { GeoElement Quadrangle  ; NumberOfPoints  4 ; }
	  { GeoElement Tetrahedron ; NumberOfPoints  4 ; }
	  { GeoElement Hexahedron  ; NumberOfPoints  6 ; }
	  { GeoElement Prism       ; NumberOfPoints  21 ; }
	  { GeoElement Line        ; NumberOfPoints  4 ; }
	}
      }
    }
  }
}

Constraint {
  { Name MVP_3D ;
    Case {
      { Region Surf_bn0 ; Type Assign ; Value 0. ; }
    }
  }
}

Group {
  Surf_a_NoGauge = Region [ {Surf_bn0, SkinDomainC} ] ;
}

Constraint {

  { Name GaugeCondition_a ; Type Assign ;
    Case {
      { Region DomainCC ; SubRegion Surf_a_NoGauge ; Value 0. ; }
    }
  }
}

FunctionSpace {

  // Magnetic vector potential a (b = curl a)
  { Name Hcurl_a_3D ; Type Form1 ;
    BasisFunction {// a = a_e * s_e
      { Name se ; NameOfCoef ae ; Function BF_Edge ;
        Support Domain ; Entity EdgesOf[ All, Not SkinDomainC] ; }
      { Name se2 ; NameOfCoef ae2 ; Function BF_Edge ;
        Support Domain ; Entity EdgesOf[ SkinDomainC] ; }
    }
    Constraint {
      { NameOfCoef ae  ; EntityType EdgesOf ; NameOfConstraint MVP_3D ; }
      { NameOfCoef ae2 ; EntityType EdgesOf ; NameOfConstraint MVP_3D ; }
      { NameOfCoef ae  ; EntityType EdgesOfTreeIn ; EntitySubType StartingOn ;
        NameOfConstraint GaugeCondition_a ; }
    }
  }

  // Electric scalar potential (3D) // not used
  { Name Hregion_u_3D ; Type Form0 ;
    BasisFunction {
      { Name sr ; NameOfCoef ur ; Function BF_GroupOfNodes ;
        Support DomainC ; Entity GroupsOfNodesOf[ Surf_Elec ] ; }
    }
    GlobalQuantity {
      { Name U ; Type AliasOf        ; NameOfCoef ur ; }
      { Name I ; Type AssociatedWith ; NameOfCoef ur ; }
    }
    Constraint {
      { NameOfCoef U ; EntityType GroupsOfNodesOf ; NameOfConstraint V_3D ; }
      { NameOfCoef I ; EntityType GroupsOfNodesOf ; NameOfConstraint I_3D ; }
    }
  }
}



Formulation {

  { Name MagStaDyn_av_js0_3D ; Type FemEquation ;
    Quantity {
      { Name a  ; Type Local ; NameOfSpace Hcurl_a_3D ; }
      //{ Name xi ; Type Local ; NameOfSpace H_xi ; }

      { Name v  ; Type Local ; NameOfSpace Hregion_u_3D ; } //Massive conductor
      { Name U  ; Type Global ; NameOfSpace Hregion_u_3D [U] ; }
      { Name I  ; Type Global ; NameOfSpace Hregion_u_3D [I] ; }
    }

    Equation {
      Galerkin { [ nu[] * Dof{d a} , {d a} ] ;
        In Domain ; Jacobian Vol ; Integration II ; }

/*       Galerkin { DtDof[ sigma[] * Dof{a} , {a} ] ;
        In DomainC ; Jacobian Vol ; Integration II ; }
      Galerkin { [ sigma[] * Dof{d v}/SymmetryFactor , {a} ] ;
        In DomainC ; Jacobian Vol ; Integration II ; }

      Galerkin { DtDof[ sigma[] * Dof{a} , {d v} ] ;
        In DomainC ; Jacobian Vol ; Integration II ; }
      Galerkin { [ sigma[] * Dof{d v}/SymmetryFactor , {d v} ] ;
        In DomainC ; Jacobian Vol ; Integration II ; }
      GlobalTerm { [ Dof{I}*SymmetryFactor, {U} ] ; In Surf_Elec ; }
 */
      Galerkin { [ -js0[], {a} ] ;
        In  DomainS ; Jacobian Vol ; Integration II ; }

    }
  }

}

Resolution {

  { Name Analysis ;
    System {
      { Name Sys ; NameOfFormulation MagStaDyn_av_js0_3D ; } // Just static
      //{ Name Sys ; NameOfFormulation MagStaDyn_av_js0_3D ; Type ComplexValue ; Frequency Freq ; }
    }
    Operation {
      CreateDir["res/"];

      InitSolution[Sys];
      Generate[Sys] ; Solve[Sys] ; SaveSolution[Sys];
      PostOperation[Get_LocalFields] ;
    //PostOperation[Get_GlobalQuantities] ;
    }
  }
}

PostProcessing {
/*   { Name TangentVector ; NameOfFormulation FindVector ;
    PostQuantity {
      { Name vdir; Value { Term { [ vDir[] ] ; In Domain ; Jacobian Vol ; } } }
    }
  } */
  { Name MagStaDyn_av_js0_3D ; NameOfFormulation MagStaDyn_av_js0_3D ;
    PostQuantity {
      { Name a ; Value { Term { [ {a} ]          ; In Domain ; Jacobian Vol ; } } }
      { Name b ; Value { Term { [ {d a} ]        ; In Domain ; Jacobian Vol ; } } }
      { Name h ; Value { Term { [ nu[] * {d a} ] ; In Domain ; Jacobian Vol ; } } }

      { Name v ; Value { Term { [ {v} ]          ; In DomainC ; Jacobian Vol ; } } }
      { Name e ; Value { Term { [ -(Dt[{a}]+{d v}) ] ; In DomainC ; Jacobian Vol ; } } }
      { Name j ; Value { Term { [ -sigma[]*(Dt[{a}]+{d v}) ] ; In DomainC ; Jacobian Vol ; } } }
      { Name js ; Value { Term { [ js0[] ]      ; In DomainS ; Jacobian Vol ; } } }

/*       { Name JouleLosses ;
        Value { Integral {
            [ SymmetryFactor * sigma[]*SquNorm[Dt[{a}]+{d v}] ] ;
            In DomainC ; Jacobian Vol ; Integration II ; }
        }
      } */

	  { Name JouleLosses ;
        Value { Integral {
            [ 0.5*SquNorm[ js0[] ]/sigma[] ] ;
            In Winding ; Jacobian Vol ; Integration II ; }
        }
      }



      { Name MagEnergy ;
        Value { Integral {
            [ 1/2 * nu[]*{d a} * {d a} ] ;
	    In Domain ; Jacobian Vol ; Integration II ; }
	}
      }

      // { Name Flux ; Value {
          // Integral { [ vDir[]/SurfCoil[]*{a} ] ;
            // In Domain  ; Jacobian Vol ; Integration II ; }
        // }
      // }

      { Name Upos ;
        Value { Integral { Type Global ;
            [ -sigma[] * (Dt[{a}] + {d v}) * BF{d v} ] ;
            In DomainC ; Jacobian Vol ; Integration II ;
          }
        }
      }


      //{ Name Inductance_from_Flux ; Value { Term { Type Global; [ $Flux * 1e3/Ipk[] ] ; In DomainDummy ; } } }
      { Name Inductance_from_Flux ; Value { Term { Type Global; [ $Flux * 1e3/Ipk[] ] ; In DomainS ; } } }
      //{ Name Inductance_from_MagEnergy ; Value { Term { Type Global; [ 2 * $MagEnergy * 1e3/(Ipk[]*Ipk[]) ] ; In DomainDummy ; } } }
      { Name Inductance_from_MagEnergy ; Value { Term { Type Global; [ 2 * $MagEnergy * 1e3/(Ipk[]*Ipk[]) ] ; In DomainS ; } } }

	  { Name Resistance ; Value { Term { Type Global; [ 2 * $JouleLosses /(Ipk[]*Ipk[]) ] ; In DomainS ; } } }
      //{ Name xi ; Value { Term { [ {xi} ] ; In Domain ; Jacobian Vol ; } } }
    }
  }
}

//========================================================================================
/*  PostOperation ShowTangentVector UsingPost TangentVector {
// added
   Print[ vdir, OnElementsOf Domain, File StrCat[Dir,"vdir",ExtGmsh] ] ;
 } */

 PostOperation Get_LocalFields UsingPost MagStaDyn_av_js0_3D {
   Print[ js, OnElementsOf DomainS, File StrCat[Dir, "js", ExtGmsh], LastTimeStepOnly ] ;
   Print[ a, OnElementsOf Domain, File StrCat[Dir, "a", ExtGmsh], LastTimeStepOnly ] ;
   // If(Flag_GaugeType==COULOMB_GAUGE)
     // Print[ xi, OnElementsOf Domain, File StrCat[Dir, "xi",ExtGmsh ], LastTimeStepOnly ] ;
   // EndIf
   Print[ b, OnElementsOf Domain, File StrCat[Dir,"b",ExtGmsh], LastTimeStepOnly ] ;
 }

  PostOperation Get_GlobalQuantities UsingPost MagStaDyn_av_js0_3D {
/*    Print[ Flux[DomainS], OnGlobal, Format TimeTable,
     File > StrCat[Dir,"Flux",ExtGnuplot], LastTimeStepOnly, StoreInVariable $Flux,
     SendToServer StrCat[po,"40Flux [Wb]"],  Color "LightYellow" ]; */

/*    Print[ Inductance_from_Flux, OnRegion Primary, Format Table, LastTimeStepOnly,
    File StrCat[Dir,"InductanceF",ExtGnuplot],
    SendToServer StrCat[po,"50Inductance from Flux [mH]"], Color "LightYellow" ]; // DomainDummy */

   Print[ MagEnergy[Domain], OnGlobal, Format TimeTable,
     File > StrCat[Dir,"ME",ExtGnuplot], LastTimeStepOnly, StoreInVariable $MagEnergy,
     SendToServer StrCat[po,"41Magnetic Energy [W]"],  Color "LightYellow" ];

   Print[ Inductance_from_MagEnergy, OnRegion Primary, Format Table, LastTimeStepOnly,
     File StrCat[Dir,"InductanceE",ExtGnuplot],
     SendToServer StrCat[po,"51Inductance from Magnetic Energy [mH]"], Color "LightYellow" ];

   Print[ JouleLosses[Primary], OnGlobal, Format TimeTable,
     File > StrCat[Dir,"WindingLoss",ExtGnuplot], LastTimeStepOnly, StoreInVariable $JouleLosses,
     SendToServer StrCat[po,"61Joule loss [W]"],  Color "LightYellow" ];

   Print[ Resistance, OnRegion Primary, Format Table, LastTimeStepOnly,
     File StrCat[Dir,"ResistanceJ",ExtGnuplot],
     SendToServer StrCat[po,"71Winding resistance [Ohm]"], Color "LightYellow" ];

}
