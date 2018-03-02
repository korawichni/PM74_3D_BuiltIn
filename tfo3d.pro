Include "tfo3d_data.geo";


// gmsh tfo3d.pro -open tfo3d_builtin.geo

Dir="res/";
ExtGmsh = ".pos";
ExtGnuplot  = ".dat";
po = "Output/Stranded/"; // For results
pcp = "Output/Circuit/Primary/"; // For results
pcs = "Output/Circuit/Secondary/"; // For results


STRANDED=0;
MASSIVE =1;

DIVJ0_NONE = 0;
DIVJ0_NONE = 0;
DIVJ0_WEAK = 1;
DIVJ0_STRONG = 2; // => to add, implementation phase

STA=0;
DYN=1;

simp = Str["Simulation parameters./"];
simlw = Str["Litz wire parameters./"];

DefineConstant[
  _winding_model = {0, Choices{0="stranded",1="massive"},
    Name StrCat[simp,"00Winding model"], Highlight "Blue"}

  _flag_litz = {_winding_model==MASSIVE ? 0 : 1, Choices{0,1}, Name StrCat[simp,"00litz wire"], ReadOnly _winding_model}	
	
  //Fill = {AreaCond/AreaCell, Name StrCat[ mgeo,"30Fill factor"], ReadOnly 1, Highlight "Blue"}
  
  _analysis_type = {_winding_model==MASSIVE ? DYN : STA, Choices{0="magnetostatics",1="magnetodynamic"},
    Name StrCat[simp,"00Choose analysis type"], Highlight "Blue", ReadOnly (_winding_model==MASSIVE)}

  _divJ_zero = { DIVJ0_STRONG,
    Choices{ DIVJ0_NONE   = "none",
             DIVJ0_WEAK   = "weak",
             DIVJ0_STRONG = "strong"},
    Name StrCat[simp,"01Constraint div j = 0"],
    Help Str["None: direct interpolation of js0[]",
      "Weak: Use scalar potential xis for weakly ensuring div j = 0.",
      "Strong: Use Hcurl source field hs with curl hs = j, for div j = 0;"],
    Highlight "Blue",  Visible (_winding_model==STRANDED)}
	
  _flag_circuit_coupling = {1, Choices{0,1}, Name StrCat[simp,"02Circuit Coupling"]}
  
];

If (_flag_litz)
  DefineConstant[
	strand_dia_pri = { strand_dia_pri_00, Name StrCat[simlw, "/05Strand diameter (primary) [m]" ], Highlight "LightGreen", ReadOnly 1, Closed !_flag_litz }
	strand_number_pri = { strand_number_pri_00, Name StrCat[simlw, "/06Number of strands in a  primary bundle" ], Highlight "LightGreen", ReadOnly 1, Visible _flag_litz } 
	strand_dia_sec = { strand_dia_sec_00, Name StrCat[simlw, "/07Strand diameter (secondary) [m]" ], Highlight "LightGreen", ReadOnly 1, Visible _flag_litz } 	
	strand_number_sec = { strand_number_sec_00, Name StrCat[simlw, "/08Number of strands in a secondary bundle" ], Highlight "LightGreen", ReadOnly 1, Visible _flag_litz }
    fill_pri = { strand_number_pri*(strand_dia_pri/(2*rp))^2, Name StrCat[simlw,"/12Primary fill factor"], Visible _flag_litz, ReadOnly 1, Highlight "LightGreen"}
	fill_sec = { strand_number_sec*(strand_dia_sec/(2*rs))^2, Name StrCat[simlw,"/13Secondary fill factor"], Visible _flag_litz, ReadOnly 1, Highlight "LightGreen"}
  ];
  
  
  
  FileRP = Sprintf("coeff/SkinProxFactors_RH_la%.2g.pro", fill_pri);//Sprintf("Coeffs_RH_la%.2g.pro", lambda_pri); "SkinProxFactors_RH_la%.2g_meeker.pro"
  FileRS = Sprintf("coeff/SkinProxFactors_RH_la%.2g.pro", fill_sec);
  Include FileRP;
  Include FileRS;
Else
  fill_pri = 1;
  fill_sec = 1;
EndIf


Group{

  Core = #CORE;
  Air  = #AIR;

  Primary    = #PRIMARY;
  Secondary0 = #SECONDARY0;
  Secondary1 = #SECONDARY1;
  Secondary   = Region[{Secondary0, Secondary1}];

  Skin_Primary    = #{SKIN_PRIMARY};
  Skin_Secondary0 = #{SKIN_SECONDARY0};
  Skin_Secondary1 = #{SKIN_SECONDARY1};

  Surf_P_In  = #{IN_PRI};
  Surf_S0_In = #{IN_SEC0};
  Surf_S1_In = #{IN_SEC1};

  Surf_P_Out  = #{OUT_PRI};
  Surf_S0_Out = #{OUT_SEC0};
  Surf_S1_Out = #{OUT_SEC1};

  // Generated automatically with Homology
  Cut_Primary    = #{1002102}; // Default: pri-1 sec0-2 sec1-3
  Cut_Secondary0 = #{1002101};
  Cut_Secondary1 = #{1002103};
 

  Surf_In  = Region[{Surf_P_In,  Surf_S0_In,  Surf_S1_In}];
  Surf_Out = Region[{Surf_P_Out, Surf_S0_Out, Surf_S1_Out}];

  Winding      = Region[{Primary, Secondary0, Secondary1}];
  Skin_Winding = Region[{Skin_Primary, Skin_Secondary0, Skin_Secondary1}];

  DomainCC = Region[ {Air, Core} ];

  If(_winding_model == STRANDED)
    DomainC = Region[{}]; // conducting domain
    SkinDomainC = Region[{}];
    Surf_Elec = Region[{}];

    If(_divJ_zero == DIVJ0_NONE || _divJ_zero == DIVJ0_WEAK)
      NbSrc_DomainB = 0;
      DomainS = Region[{Winding}];
      SkinDomainS = Region[{Skin_Winding}];

      DomainB     = Region[{}];// Using js0[] directly
      SkinDomainB = Region[{}];
    EndIf
    If(_divJ_zero == DIVJ0_STRONG)
      NbSrc_DomainB = 3;
      DomainB   = Region[{Winding}];// Using {d hs}, precomputation of hs
      SkinDomainB = Region[{Skin_Winding}];

      DomainS   = Region[{}];
      SkinDomainS = Region[{}];
    EndIf
  Else
    DomainC = Region[{Winding}];
    SkinDomainC = Region[{Skin_Winding}];

    Surf_Elec = Region[{Surf_In}];
    DomainS   = Region[{}];
    SkinDomainS  = Region[{}];

    NbSrc_DomainB = 0;
    DomainB = Region[{}];
    SkinDomainB = Region[{}];
  EndIf
  DomainCC += Region[{DomainS, DomainB}];

  Surf_bn0 = #SURF_AIROUT;
  Surf_FixedMVP = Region[{Surf_bn0, Surf_In, Surf_Out}];


  // Source: stranded coils
  SurfGh0 = Region[{}];

  DomainB~{1}    = Region[{Primary}];
  SurfCutB~{1}   = Region[{Cut_Primary}] ;
  SkinDomainB~{1}= Region[{Skin_Primary}];

  DomainB~{2}    = Region[{Secondary0}];
  SurfCutB~{2}   = Region[{Cut_Secondary0}] ;
  SkinDomainB~{2}= Region[{Skin_Secondary0}];

  DomainB~{3}    = Region[{Secondary1}];
  SurfCutB~{3}   = Region[{Cut_Secondary1}] ;
  SkinDomainB~{3}= Region[{Skin_Secondary1}];

  For k In { 1:NbSrc_DomainB }
    SkinDomainB_tot~{k} = Region[ {SkinDomainB~{k}, SurfGh0} ] ;
  EndFor

  Domain = Region[ {DomainC, DomainCC} ];
  DomainDummy = Region[ {123456789} ];
}

Function{
  // RMS current
  Irms1 = 0;
  Irms2 = 10;

  // To be adapted
  If (_flag_litz)
    Nw_pri  = strand_number_pri;//1.;
    Nw_sec0 = strand_number_sec;
    Nw_sec1 = strand_number_sec;
  Else
	Nw_pri  = 1;//1.;
    Nw_sec0 = 1;
    Nw_sec1 = 1;  
  EndIf
	
  DefineConstant[
    Freq = { 50., Min 0, Max 1e3, Step 1,
      Name StrCat[simp,"10source/00frequency [Hz]"], Highlight "AliceBlue"},
    Irms_pri = { Irms1, Min 1, Max 4*Irms1, Step 2,
      Name StrCat[simp,"10source/01primary current (rms) [A]"], Highlight "AliceBlue"}
    Irms_sec = { Irms2, Min 1, Max 4*Irms2, Step 2,
      Name StrCat[simp,"10source/02secondary current (rms) [A]"], Highlight "AliceBlue"}

    sigma_coil = { sigma_cu,
      Name  StrCat[simp,"10source/10conductivity"], Units "S/m", Highlight "AliceBlue"},
    mur_fe = { 2000, Min 100, Max 2000, Step 100,
      Name  StrCat[simp,"11core relative permeability"], Highlight "AliceBlue"}
    //N27 initial permeability = 2000

    SymmetryFactor = 1
  ];

  NbWires[Primary]    = Nw_pri;
  NbWires[Secondary0] = Nw_sec0;
  NbWires[Secondary1] = Nw_sec1;

  // Peak currents for AC case; For DC case, using the rms signal or DC signal yielding the same average power
  IA_pri  = Irms_pri * Sqrt[2];
  IA_sec0 = Irms_sec * Sqrt[2];
  IA_sec1 = Irms_sec * Sqrt[2];

  IA[#{Primary, Surf_P_In}]     = IA_pri ; 
  IA[#{Secondary0, Surf_S0_In}] = IA_sec0;
  IA[#{Secondary1, Surf_S1_In}] = IA_sec1;
  

  IA[#{Air,Core}] = IA_pri*0;

  //FactorAC = (_analysis_type==STA)? 1 : 2; // In STA mode, we inject the peak signal, in DC, P = Ipk*Vpk, in AC, Sav = 0.5*Ipk*Vpk
  Ap = (interwire_pri+2*(rp+thick_insul))/4 ;
  As = (interwire_sec+2*(rs+thick_insul))/4 ;
/*
  Np_=Np-0.25*0;
  vDir[#{Primary}] =
  (Y[] >= yp0-rp && Z[] >=0 && X[]>0) ? Vector [0.,0.,-1.] :
  (Y[] <= yp0+2*rp -(interwire_pri+2*(rp+thick_insul))*Np_ && Z[]>=0 && X[]>0) ? Vector [1.,0.,0.]:
  Unit[ Vector[Z[], -Ap ,-X[]] ] ;

  Ns_=Ns-0.25;
  vDir[#{Secondary0,Secondary1}] =
  ( (Y[] >= ys0-rs && Z[] >=0 && X[]>0) ? Vector [0.,0.,-1.]:
    (Y[] <= ys0+2*rs-(interwire_sec+2*(rs+thick_insul))*Ns_ && Z[]>=0 && X[] >=0) ? Vector [1.,0.,0.]:
    Unit [ Vector [ Z[], -As ,-X[] ] ] ); */
  //================================================
  /* Parameters w and a depending on how to construct the winding
         x(t) = R*cos(wt), w = Pi/2
         y(t) = H - at,    a = pitch/4
	     z(t) = -R*sin(wt)
  The tangent vector of r(t) = {x(t),y(t),z(t)} is r'(t) = {w*z(t), -a, -w*x(t)} being normalized 
  *///===============================================

  
  eps = 1e-3*1;
  _tmp_omega = Pi/2;
  pitch_pri = (interwire_pri+2*(rp+thick_insul));
  pitch_sec = (interwire_sec+2*(rs+thick_insul));
  _tmp_pri = yp0+rp-pitch_pri*(Np-0.25) + eps;  
  _tmp_sec = ys0+rs-pitch_sec*(Ns-0.25) + eps; // this value is -14.195 mm in raw value.

  vDir[Primary] = ( (Y[] >= yp0-rp && Z[] >=0 && X[] >=0 ) ? Vector [0,0,-1]: 
                  ( Y[]<= _tmp_pri  && X[] >=0 && Z[] >=0) ? Vector [1,0,0]: 
                    Unit [ Vector [ _tmp_omega*Z[], -Ap ,-_tmp_omega*X[] ] ] );
  vDir[#{Secondary0,Secondary1}] = ( (Y[] >= ys0-rs && Z[] >=0 && X[] >=0 ) ? Vector [0,0,-1]:
                    ( Y[]<=  _tmp_sec && X[] >=0 && Z[] >=0) ? Vector [1,0,0]:
                      Unit [ Vector [ _tmp_omega*Z[], -As ,-_tmp_omega*X[] ] ] );
	
	
  SurfCoil[Primary]    = SurfaceArea[]{IN_PRI};
  SurfCoil[Secondary0] = SurfaceArea[]{IN_SEC0};
  SurfCoil[Secondary1] = SurfaceArea[]{IN_SEC1};

  js1A[] = NbWires[]/SurfCoil[]*vDir[];
  
  js0[]  = _flag_litz? IA[]/NbWires[]*js1A[] : IA[]*js1A[];

  // Material properties
  mu0 = 4.e-7 * Pi ;
  nu0 = 1./mu0;  
  nu[ Region[{Air}]] = nu0;
  nu[Core] = 1/(mur_fe*mu0) ;

  sigma[Winding] = sigma_coil ;
  rho[] = 1/sigma[] ;

  Rdc_pri =  (Hypot[2*Pi*xp0,pitch_pri]*(Np-0.25) + 2*z0)/sigma_coil/(fill_pri*Pi*(rp^2));
  Rdc_sec0 =  (Hypot[2*Pi*xs0,pitch_sec]*(Ns-0.25) + 2*z0)/sigma_coil/(fill_sec*Pi*(rs^2));
  Rdc_sec1 =  (Hypot[2*Pi*xs1,pitch_sec]*(Ns-0.25) + 2*z0)/sigma_coil/(fill_sec*Pi*(rs^2));
  Rdc_sec = Rdc_sec0 + Rdc_sec1;
  Printf("Rdc_pri=",Rdc_pri); //0.00180492 Ohm
  Printf("Rdc_sec=",Rdc_sec); //0.0102548 Ohm
  // Homogenization coefficients: round conductor & hexagonal packing 
  If (_flag_litz)
    DefineConstant[
      delta = {Sqrt[1/(Pi*Freq*mu0*sigma_coil)] , Name StrCat[simlw,"/20Skin depth (m)"], Visible (_analysis_type==DYN) ,Highlight "LightGreen", ReadOnly 1 }	
      X_pri = { strand_dia_pri*1/2/delta, Name StrCat[simlw,"/21Reduced frequncy ratio, Xpri"], Visible (_analysis_type==DYN), ReadOnly 1, Highlight "LightGreen"}
      X_sec = { strand_dia_sec*1/2/delta, Name StrCat[simlw,"/22Reduced frequncy ratio, Xsec"], Visible (_analysis_type==DYN), ReadOnly 1, Highlight "LightGreen"}
      _Rdc_pri = { Rdc_pri, Name StrCat[simlw,"/23DC resistance (pri) [Ohm]"], Visible (_analysis_type==DYN), ReadOnly 1, Highlight "LightGreen"}	  
      _Rdc_sec = { Rdc_sec, Name StrCat[simlw,"/23DC resistance (sec) [Ohm]"], Visible (_analysis_type==DYN), ReadOnly 1, Highlight "LightGreen"}	  
    ];
    
	// Frequency domain
    skin_rhor_pri[] = InterpolationLinear[$1]{List[skin_rhor_Circ_pri]} ;
    skin_rhoi_pri[] = InterpolationLinear[$1]{List[skin_rhor_Circ_pri]} ;
    prox_nur_pri[]  = InterpolationLinear[$1]{List[prox_nur_Circ_pri]} ;
    prox_nui_pri[]  = InterpolationLinear[$1]{List[prox_nui_Circ_pri]} ;

    skin_rhor_sec[] = InterpolationLinear[$1]{List[skin_rhor_Circ_sec]} ;
    skin_rhoi_sec[] = InterpolationLinear[$1]{List[skin_rhor_Circ_sec]} ;
    prox_nur_sec[]  = InterpolationLinear[$1]{List[prox_nur_Circ_sec]} ;
    prox_nui_sec[]  = InterpolationLinear[$1]{List[prox_nui_Circ_sec]} ;

	fill[Primary] = fill_pri;
	fill[Secondary] = fill_sec;
	
    nu [Primary] = nu0*Complex[ prox_nur_pri[X_pri], prox_nui_pri[X_pri]*fill_pri*X_pri^2/2];//* Complex[qb_pri,pb_pri*lambda_pri*x_pri^2/2];
    nu [Secondary] = nu0*Complex[ prox_nur_sec[X_sec], prox_nui_sec[X_sec]*fill_sec*X_sec^2/2];//*
	
	// zskin_pri[Primary] = Rdc_pri*Complex[ skin_rhor[X_pri],skin_rhoi[X_pri]*X_pri^2/4/fill_pri];
	// zskin_sec[Secondary] = Rdc_sec*Complex[ skin_rhor[X_sec],skin_rhoi[X_sec]*X_sec^2/4/fill_sec];
    AreaCell[Primary] = strand_number_pri;///strand_number_pri;//Pi*rp^2;   //-1/(strand_number_pri*(Pi*(strand_dia_pri/2)^2)/fill_pri);//1/SurfaceArea[];//NoOfPri/SurfaceArea[]; // number of turns per surface area (only + sign)
    AreaCell[Secondary] =  strand_number_sec;//Pi*rs^2;  //-1/(strand_number_sec*(Pi*(strand_dia_sec/2)^2)/fill_sec);//1/SurfaceArea[];//NoOfSec/SurfaceArea[]; // number of turns per surface area (only + sign)
 
  Else 
   fill[#{Primary,Secondary}] = 1;
   nu[ Region[{Winding}]] = nu0;
   AreaCell[#{Primary,Secondary}] = 1;
  EndIf
 
}


Jacobian {
  { Name Vol ;
    Case {
      { Region All ; Jacobian Vol ; }
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
      {	Type Gauss ;
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

// ---------
//  Circuit
// ---------
If (_flag_circuit_coupling)
	
  Group {
    VI_source_pri = # 2000001 ;
    Rs_pri  = # 2000002 ;
	Rz_pri = # 2000003 ;
	
    VI_source_sec = # 2000004 ;
    Rs_sec  = # 2000005 ;
	Rz_sec = # 2000006 ;
	
    //Resistance_Cir  = Region[ {Rs_pri,Rz_pri,Rs_sec,Rz_pri} ] ;
    Resistance_Cir  = Region[ {Rs_pri,Rs_sec,Rz_pri,Rz_sec} ] ;
    Inductance_Cir  = Region[ {} ] ;

    Capacitance1_Cir = Region[ {} ] ;
    Capacitance2_Cir = Region[ {} ] ;
    Capacitance_Cir = Region[ {Capacitance1_Cir, Capacitance2_Cir} ] ;

    Diode_Cir  = Region[ {} ] ;

    SourceV_Cir = Region[ { } ] ;
    //SourceI_Cir = Region[ {VI_source_pri, VI_source_sec} ] ;
    SourceI_Cir = Region[ {VI_source_pri, VI_source_sec} ] ;

    DomainZ_Cir = Region[ {Resistance_Cir, Inductance_Cir, Capacitance_Cir,
                           Diode_Cir} ] ;
    DomainSource_Cir = Region[ {SourceV_Cir, SourceI_Cir} ] ;

    DomainZt_Cir = Region[ {DomainZ_Cir, DomainSource_Cir} ] ;
  }
  
  Function {
	DefineFunction[ Resistance ];
	// If (!Flag_Ic)
	  // Resistance[Rc] = Rload ;
	// EndIf
	// If (Flag_Ic)
	Resistance[#{Rs_pri,Rs_sec}] = 1e50;
    If (_flag_litz)
  	  Resistance[Rz_pri] = Rdc_pri*Complex[ skin_rhor_pri[X_pri],skin_rhoi_pri[X_pri]*X_pri^2/4/fill_pri];
  	  Resistance[Rz_sec] = Rdc_sec*Complex[ skin_rhor_sec[X_sec],skin_rhoi_sec[X_sec]*X_sec^2/4/fill_sec]; 
	Else
	  Resistance[#{Rz_pri,Rz_sec}] = 0;
	EndIf
	//EndIf
  }
  
  Constraint {
    { Name Current_Cir ;
    Case {
	If(_analysis_type == STA) // Static and imposed current
		{ Region VI_source_pri ; Value Irms_pri ;}
		{ Region VI_source_sec ; Value Irms_sec ;}
	EndIf
	If(_analysis_type == DYN) // sinusoidal voltage source, either TL or FD   
      { Region VI_source_pri ; Value IA_pri ; TimeFunction F_Cos_wt_p[]{2*Pi*Freq,0} ; } 
	  { Region VI_source_sec ; Value IA_sec0 ; TimeFunction F_Cos_wt_p[]{2*Pi*Freq,0} ; } 
	EndIf
      }
    }
	
    { Name Voltage_Cir ; Case { } }	

  // { Name Voltage_Cir ;
    // Case {
	// If(_analysis_type == STA) // Static and imposed voltage
		// { Region VI_Source ; Value Vc ;}
	// EndIf
	// If(Flag_TL && !Flag_SinusVI && Flag_Vc /*Flag_TL_Vstep*/) // TL and voltage step
        // { Region VI_Source ; Value Vquasi_sq_TL ; TimeFunction Quasi_sq[] ;  }         
    // EndIf
	// If(_analysis_type == DYN) // sinusoidal voltage source, either TL or FD   
      // { Region VI_Source ; Value Vpri_pk ; TimeFunction F_Cos_wt_p[]{2*Pi*Freq,Vc_ph*Pi/180} ; } 
	// EndIf
    // }
  // }

    { Name ElectricalCircuit ; Type Network ;

//      1 _________ -- VI_source_pri --> _________ 3
//      |                                          |
//      |_________ Rz_pri ______2____ Primary _____|
//		|							               |
//		|___________________Rs_pri_________________|

	  Case Circuit1 {
		If (Np > 0)
		  { Region Primary ; Branch {2, 3} ; }
		EndIf 
		
		{ Region Rz_pri ; Branch {1, 2} ; }
		{ Region VI_source_pri ; Branch {1, 3} ; }
		{ Region Rs_pri; Branch {1, 3} ; }
	  }
	
//      1 ______________ -- VI_source_sec --> _____________ 4
//      |                                                   |
//      |___ Rz_sec ___2___ Secondary0 __3__ Secondary1_____|
//		|							                        |
//		|______________________Rs_sec_______________________|
	
      Case Circuit2 {
		If (Ns > 0)
		  { Region Secondary0 ; Branch {2, 3} ; }
		  { Region Secondary1 ; Branch {3, 4} ; }
		EndIf
		
		{ Region Rz_sec ; Branch {1, 2} ; }		
		{ Region VI_source_sec ; Branch {1, 4} ; }
		{ Region Rs_sec; Branch {1, 4} ; } // Need this to make it well-conditioned
	  }	 
   } 
  } // End of Constraint
  
  // UZ and IZ for impedances
  FunctionSpace {
    { Name Hregion_Z ; Type Scalar ;
      BasisFunction {
        { Name sr ; NameOfCoef ir ; Function BF_Region ;
          Support DomainZt_Cir ; Entity DomainZt_Cir ; }
      }
      GlobalQuantity {
        { Name Iz ; Type AliasOf        ; NameOfCoef ir ; }
        { Name Uz ; Type AssociatedWith ; NameOfCoef ir ; }
      }
      Constraint {
        { NameOfCoef Uz ;
          EntityType Region ; NameOfConstraint Voltage_Cir ; }
        { NameOfCoef Iz ;
          EntityType Region ; NameOfConstraint Current_Cir ; }
      }
    }
  
  }
  
EndIf // EndIf of CircuitCoupling  

Constraint {
  { Name MVP_3D ;
    Case {
      { Region Surf_bn0 ; Type Assign ; Value 0. ; }
      { Region Surf_In  ; Type Assign ; Value 0. ; }
      { Region Surf_Out ; Type Assign ; Value 0. ; }
    }
  }

  { Name V_3D ;
    Case {
	  If (!_flag_circuit_coupling)
        If(_winding_model==STRANDED)
		  If (IA_pri == 0)
		    { Region Primary ; Type Assign ; Value 0. ; }
		  EndIf
		  If (IA_sec0 == 0)
			{ Region Secondary0 ; Type Assign ; Value 0. ; }
            { Region Secondary1 ; Type Assign ; Value 0. ; }
		  EndIf
        EndIf
	  EndIf
    }
  }

  { Name I_3D ;
    Case {
	  If (!_flag_circuit_coupling)
        If(_winding_model==MASSIVE)
          { Region Surf_P_In  ; Type Assign ; Value IA[] ; TimeFunction 1.; }
          { Region Surf_S0_In ; Type Assign ; Value IA[] ; TimeFunction 1.; }
          { Region Surf_S1_In ; Type Assign ; Value IA[] ; TimeFunction 1.; }
        EndIf
        If(_winding_model==STRANDED)
	      If (IA_pri != 0)
            { Region Primary ; Type Assign ; Value /*-IA[]*/ Irms_pri; TimeFunction 1.; }
		  EndIf
	      If (IA_sec0 != 0)
            { Region Secondary0 ; Type Assign ; Value Irms_sec ; TimeFunction 1.; } 
            { Region Secondary1 ; Type Assign ; Value Irms_sec ; TimeFunction 1.; }
		  EndIf		  
        EndIf
	  EndIf
    }
  }

  // Constraints for computing source magnetic field hs
  { Name I_Unit ;
    Case {
      { Region SurfCutB~{1}; Value Nw_pri ; }
      { Region SurfCutB~{2}; Value Nw_sec0 ; }
      { Region SurfCutB~{3}; Value -Nw_sec1 ; }
    }
  }

  For k In { 1:NbSrc_DomainB }
    { Name GaugeCondition_hs~{k} ; Type Assign ;
      Case {
        { Region DomainB~{k} ; SubRegion SkinDomainB_tot~{k} ; Value 0. ; }
      }
    }
  EndFor

  { Name MagneticField ;
    Case {
    }
  }

}

//------------------------------------------------------------------------------------
// -- For computing: Source magnetic field -- stranded coils (bobine)
//------------------------------------------------------------------------------------
For k In {1:NbSrc_DomainB}
  FunctionSpace{
    // Stranded conductors: magnetic field hs (generalized source field) -- static
    { Name Hcurl_hs~{k} ; Type Form1 ;
      BasisFunction {
        { Name se ; NameOfCoef he ; Function BF_Edge ;
          Support DomainB~{k} ; Entity EdgesOf[ All, Not SkinDomainB_tot~{k} ] ; }
        // { Name sc ; NameOfCoef Ic ; Function BF_GradGroupOfNodes ;
        //  Support ElementsOf[DomainB~{k}, OnOneSideOf SurfCutB~{k} ] ;
        //  Entity GroupsOfNodesOf[ SurfCutB~{k} ] ; }
        { Name sc ; NameOfCoef Icc ; Function BF_GroupOfEdges ;
          Support DomainB~{k} ; Entity GroupsOfEdgesOf[ SurfCutB~{k} ] ; }
      }

      Constraint {
        { NameOfCoef he  ; EntityType EdgesOf ; NameOfConstraint MagneticField ; }
        { NameOfCoef he ;  // Gauge condition
          EntityType EdgesOfTreeIn ; EntitySubType StartingOn ;
          NameOfConstraint GaugeCondition_hs~{k} ; }
        { NameOfCoef Icc ; EntityType GroupsOfNodesOf ; NameOfConstraint I_Unit ; }
      }
    }
  }

  Formulation {
    { Name MagSta_hs~{k}; Type FemEquation;
      Quantity {
        { Name hs; Type Local; NameOfSpace Hcurl_hs~{k}; }
      }
      Equation {
        Galerkin { [ Dof{d hs} , {d hs} ] ;
          In DomainB~{k}; Integration II; Jacobian Vol;  }
        Galerkin { [   -js1A[] , {d hs} ] ;
          In DomainB~{k}; Integration II; Jacobian Vol;  }
      }
    }
  }

  Resolution {
    { Name MagSta_hs~{k} ;
      System {
        { Name H ; NameOfFormulation MagSta_hs~{k} ; }
      }
      Operation {
        Generate[H] ; Solve[H] ;
        // PostOperation[Map_hs~{k}];
        // SaveSolution[H] ;
      }
    }
  }

  PostProcessing {
    { Name MagSta_hs~{k}; NameOfFormulation MagSta_hs~{k}; NameOfSystem H;
      Quantity {
        { Name hs;  Value{ Local{ [ {hs} ]   ; In DomainB~{k}; Jacobian Vol; } } }
        { Name js;  Value{ Local{ [ {d hs} ] ; In DomainB~{k}; Jacobian Vol; } } }
        { Name js1A; Value{ Local{ [ js1A[] ]  ; In DomainB~{k}; Jacobian Vol; } } }
      }
    }
  }

  PostOperation {
    { Name Map_hs~{k} ; NameOfPostProcessing MagSta_hs~{k} ;
      Operation {
        Print[ hs,   OnElementsOf DomainB~{k}, File Sprintf("hs_%g.pos",k) ] ;
        Print[ js,   OnElementsOf DomainB~{k}, File Sprintf("js_%g.pos",k) ] ;
        Print[ js1A, OnElementsOf DomainB~{k}, File Sprintf("js1A_%g.pos",k)] ;
      }
    }
  }
EndFor
//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

Group {
  Surf_a_NoGauge = Region [ {Surf_FixedMVP, SkinDomainC} ] ;
}

Constraint {

  { Name GaugeCondition_a ; Type Assign ;
    Case {
      { Region DomainCC ; SubRegion Surf_a_NoGauge ; Value 0. ; }
    }
  }

  { Name xi_fixed ; Type Assign ;
    Case {
      { Region Surf_FixedMVP ; Value 0. ; }
      { Region SkinDomainC   ; Value 0. ; }
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

  // Electric scalar potential (3D), used if massive conductor with applied I or V
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


  // Source field in AV formulation for stranded conductors
  { Name HsSpace ; Type Form1 ;
    BasisFunction {
      { Name sr ; NameOfCoef ir ;  // Global Basis Function
        Function BF_Global { Quantity hs ;
          Formulation MagSta_hs {NbSrc_DomainB} ;
          Group DomainB ; Resolution MagSta_hs {NbSrc_DomainB} ; } ;
        Support Domain ; Entity Global [DomainB] ; }
    }
    GlobalQuantity {
      { Name Ib ; Type AliasOf        ; NameOfCoef ir ; }
      { Name Ub ; Type AssociatedWith ; NameOfCoef ir ; }
    }
    Constraint {
      { NameOfCoef Ib ; EntityType Global ; NameOfConstraint I_3D ; }
      { NameOfCoef Ub ; EntityType Global ; NameOfConstraint V_3D ; }
    }
  }


  // correcting source interpolation js0[] so that (weakly) div j = 0
  { Name H_xi_divj0 ; Type Form0 ;
    BasisFunction {
      { Name sn ; NameOfCoef an ; Function BF_Node ;
        Support Region[{DomainS, SkinDomainS}] ; Entity NodesOf[ All ] ; }
    }
    Constraint {
      { NameOfCoef an ; EntityType NodesOf ; NameOfConstraint xi_fixed ; }
    }
  }

}



Formulation {

  { Name DivJ0 ; Type FemEquation ;
    Quantity {
      { Name xis; Type Local ; NameOfSpace H_xi_divj0 ; }
    }
    Equation {
      Galerkin { [    js0[] , {d xis} ] ;
        In Domain ; Jacobian Vol ; Integration II ; }
      Galerkin { [ -Dof{d xis} , {d xis} ] ;
        In Domain ; Jacobian Vol ; Integration II ; }
    }
  }

  { Name MagStaDyn_av_js0_3D ; Type FemEquation ;
    Quantity {
      { Name a  ; Type Local ; NameOfSpace Hcurl_a_3D ; }

      // Massive conductors (including coils)
      { Name v ; Type Local  ; NameOfSpace Hregion_u_3D    ; }
      { Name U ; Type Global ; NameOfSpace Hregion_u_3D[U] ; }
      { Name I ; Type Global ; NameOfSpace Hregion_u_3D[I] ; }

      // Stranded coils
      { Name hs ; Type Local  ; NameOfSpace HsSpace ; }
      { Name Ib ; Type Global ; NameOfSpace HsSpace[Ib] ; }
      { Name Ub ; Type Global ; NameOfSpace HsSpace[Ub] ; }

      { Name xis ; Type Local ; NameOfSpace H_xi_divj0 ; } // div j=0
	  
	  If (_flag_circuit_coupling)
        { Name Uz ; Type Global ; NameOfSpace Hregion_Z [Uz] ; }  //For the lumped components and voltage
        { Name Iz ; Type Global ; NameOfSpace Hregion_Z [Iz] ; }  //For the lumped components and current 	  
      EndIf
    }

    Equation {
      Galerkin { [ nu[] * Dof{d a} , {d a} ] ;
        In Domain ; Jacobian Vol ; Integration II ; }

      Galerkin { DtDof[ sigma[] * Dof{a} , {a} ] ;
        In DomainC ; Jacobian Vol ; Integration II ; }
      Galerkin { [ sigma[] * Dof{d v}/SymmetryFactor , {a} ] ;
        In DomainC ; Jacobian Vol ; Integration II ; }

      Galerkin { DtDof[ sigma[] * Dof{a} , {d v} ] ;
        In DomainC ; Jacobian Vol ; Integration II ; }
      Galerkin { [ sigma[] * Dof{d v}/SymmetryFactor , {d v} ] ;
        In DomainC ; Jacobian Vol ; Integration II ; }
      GlobalTerm { [ Dof{I}*SymmetryFactor, {U} ] ; In Surf_Elec ; }
	  
	  If(!_flag_circuit_coupling)
		Galerkin { [ -js0[], {a} ] ;
          In DomainS ; Jacobian Vol ; Integration II ; }
	  EndIf
	  
      If(_divJ_zero == DIVJ0_WEAK)
        Galerkin { [ {d xis}, {a} ] ;
          In Domain ; Jacobian Vol ; Integration II ; }
      EndIf
      If(_divJ_zero == DIVJ0_STRONG)
        // Stranded coil
        Galerkin { [ -1/AreaCell[]*Dof{d hs} , {a} ];
          In DomainB; Jacobian Vol ; Integration II ; }
        Galerkin { DtDof [ 1/AreaCell[]* Dof{a} , {d hs} ];
          In DomainB; Jacobian Vol ; Integration II ; }
		  
		If (!_flag_litz || !_flag_circuit_coupling) // Impose Joules'law on the material, which is required in the simple stranded model
          Galerkin { [ 1/(fill[]*sigma[]) * 1/(AreaCell[]*AreaCell[]) * Dof{d hs} , {d hs}];
            In DomainB; Jacobian Vol ; Integration II ; }
		EndIf
		
        GlobalTerm { [ Dof{Ub}/SymmetryFactor , {Ib} ] ; In DomainB ; }
      EndIf
    
	  If(_flag_circuit_coupling)
	    GlobalTerm { NeverDt[ Dof{Uz}  , {Iz} ] ; In Resistance_Cir ; }
	    GlobalTerm { NeverDt[ Resistance[] * Dof{Iz} , {Iz} ] ; In Resistance_Cir ; } // Zskin considered here

	    // GlobalTerm { [ Dof{Uz}                      , {Iz} ] ; In Inductance_Cir ; }
	    // GlobalTerm { DtDof [ Inductance[] * Dof{Iz} , {Iz} ] ; In Inductance_Cir ; }
		GlobalTerm { [ 0. * Dof{Iz} , {Iz} ] ; In DomainZt_Cir ; }
		GlobalTerm { [ 0. * Dof{Uz} , {Iz} ] ; In DomainZt_Cir ; }
	  
		// imposing the Kirchhoff nodal equations (sum of currents in each node = 0)
		// and the Kirchhoff loop equations (sum of voltages in each loop = 0)
		GlobalEquation {
		  Type Network ; NameOfConstraint ElectricalCircuit ;
		  { Node {Ib};  Loop {Ub};  Equation {Ub};  In DomainB; } // for coils domain
		  //{ Node {I};  Loop {U};  Equation {I};  In DomainC_Mag ; } // for massive conductors domain
		  { Node {Iz}; Loop {Uz}; Equation {Uz}; In DomainZt_Cir ; }
	}
      EndIf
		
	}
  }

}

Resolution {

  { Name Analysis ; // not completely general: to adapt
    System {
      If( _analysis_type==STA)
        If(_divJ_zero == DIVJ0_WEAK)
          { Name Sys_DivJ0 ; NameOfFormulation DivJ0 ; }
        EndIf
        { Name Sys ; NameOfFormulation MagStaDyn_av_js0_3D ; } // static
      Else
        If(_divJ_zero == DIVJ0_WEAK)
          { Name Sys_DivJ0 ; NameOfFormulation DivJ0 ; Type ComplexValue ; Frequency Freq ;}
        EndIf
        { Name Sys ; NameOfFormulation MagStaDyn_av_js0_3D ; Type ComplexValue ; Frequency Freq ; }
      EndIf
     }
    Operation {
      CreateDir["res/"];

      If(_divJ_zero == DIVJ0_WEAK)
        Generate[Sys_DivJ0] ; Solve[Sys_DivJ0] ; SaveSolution[Sys_DivJ0];
      EndIf
      InitSolution[Sys];
      Generate[Sys] ; Solve[Sys] ; SaveSolution[Sys];
      PostOperation[Get_LocalFields] ;
      PostOperation[Get_GlobalQuantities] ;
    }
  }
}

PostProcessing {

  { Name MagStaDyn_av_js0_3D ; NameOfFormulation MagStaDyn_av_js0_3D ;
    PostQuantity {
      { Name a ; Value { Term { [ {a} ]        ; In Domain ; Jacobian Vol ; } } }
      { Name b ; Value { Term { [ {d a} ]      ; In Domain ; Jacobian Vol ; } } }
      { Name h ; Value { Term { [ nu[{d a}]*{d a} ] ; In Domain ; Jacobian Vol ; } } }

      { Name v ; Value { Term { [ {v} ]              ; In DomainC ; Jacobian Vol ; } } }
      { Name e ; Value { Term { [ -(Dt[{a}]+{d v}) ] ; In DomainC ; Jacobian Vol ; } } }
      { Name j ; Value { Term { [ -sigma[]*(Dt[{a}]+{d v}) ] ; In DomainC ; Jacobian Vol ; } } }

      { Name js0 ; Value { Term { [ js0[] ]            ; In DomainS ; Jacobian Vol ; } } }
      { Name js  ; Value { Term { [ {d hs}/AreaCell[] ]           ; In Domain ; Jacobian Vol ; } } }
      { Name hs  ; Value { Term { [ {hs} ]             ; In Domain ; Jacobian Vol ; } } }

      { Name xis ; Value { Term { [ {xis} ]   ; In Domain ; Jacobian Vol ; } } }
      { Name dxis; Value { Term { [ {d xis} ] ; In Domain ; Jacobian Vol ; } } }
      { Name js0_dxis; Value { Term { [ js0[]-{d xis} ] ; In Domain ; Jacobian Vol ; } } }

		
      { Name JouleLosses ;
        Value {
		If (_analysis_type == MASSIVE)
          Integral {[ SymmetryFactor * sigma[]*SquNorm[Dt[{a}]+{d v}] ] ;
            In DomainC ; Jacobian Vol ; Integration II ; }
		EndIf
		If (_analysis_type == STRANDED)
			
		  If (_divJ_zero == DIVJ0_NONE)
            Integral {[ SymmetryFactor * SquNorm[ js0[] ]/sigma[] ] ;
            In Domain ; Jacobian Vol ; Integration II ; }
	      EndIf
		  
		  If (_divJ_zero == DIVJ0_WEAK)
            Integral {[ SymmetryFactor * SquNorm[ (js0[]-{d xis}) ]/sigma[] ] ;
            In Domain ; Jacobian Vol ; Integration II ; }
	      EndIf
		  
		  If (_divJ_zero == DIVJ0_STRONG)
			If (!_flag_litz)
              Integral {[ SymmetryFactor * SquNorm[ {d hs} ]/sigma[] ] ;
              In Domain ; Jacobian Vol ; Integration II ; }
			Else
              Integral {[ SymmetryFactor * SquNorm[ 1/AreaCell[] * {d hs} ]/(fill[]*sigma[]) ] ;
              In Domain ; Jacobian Vol ; Integration II ; }
			EndIf
	      EndIf
		  
		EndIf
        }
      }

      { Name MagEnergy ;
        Value { Integral {
            [ SymmetryFactor * 1/2 * nu[{d a}] * {d a} * {d a} ] ;
	    In Domain ; Jacobian Vol ; Integration II ; }
	    }
      }

      { Name Flux ; Value {
          Integral { [ SymmetryFactor*vDir[]*NbWires[]/SurfCoil[]*{a} ] ;
            In DomainS  ; Jacobian Vol ; Integration II ; }
        }
      }

      { Name Upos ;
        Value { Integral { Type Global ;
            [ -sigma[] * (Dt[{a}] + {d v}) * BF{d v} ] ;
            In DomainC ; Jacobian Vol ; Integration II ; }
        }
      }
	  

	  
  If (_flag_circuit_coupling)
	  { Name U ; Value {
          Term { [ {U} ]   ; In DomainC ; }
          Term { [ {Ub} ]   ; In DomainB ; }
          Term { [ {Uz} ]  ; In DomainZt_Cir ; }
        } 
	  }
      { Name I ; Value {
          Term { [ {I} ]   ; In DomainC ; }
          Term { [ {Ib} ]   ; In DomainB ; }
          Term { [ {Iz} ]  ; In DomainZt_Cir ; }
        }
	  }
	  
	  { Name Rac ;
      Value {
	  Term { Type Global; [ Re[{Uz}/{Iz}] ] ; In DomainZt_Cir ;  } } }	  
	  
     If (_analysis_type == DYN)
		
	  { Name Upk ; Value {
          Term { [ Norm[{U}] ]   ; In DomainC ; }
          Term { [ Norm[{Ub}] ]   ; In DomainB ; }
          Term { [ Norm[{Uz}] ]  ; In DomainZt_Cir ; }
        } 
	  }
	  { Name Uph ; Value {
          Term { [ 180/Pi*Atan2[Im[{U}],Re[{U}]] ]   ; In DomainC ; }
          Term { [ 180/Pi*Atan2[Im[{Ub}],Re[{Ub}]] ]   ; In DomainB ; }
          Term { [ 180/Pi*Atan2[Im[{Uz}],Re[{Uz}]] ]  ; In DomainZt_Cir ; }
        } 
	  }
		
	  { Name Ipk ; Value {
          Term { [ Norm[{I}] ]   ; In DomainC ; }
          Term { [ Norm[{Ib}] ]   ; In DomainB ; }
          Term { [ Norm[{Iz}] ]  ; In DomainZt_Cir ; }
        } 
	  }
	  { Name Iph ; Value {
          Term { [ 180/Pi*Atan2[Im[{I}],Re[{I}]] ]   ; In DomainC ; }
          Term { [ 180/Pi*Atan2[Im[{Ib}],Re[{Ib}]] ]   ; In DomainB ; }
          Term { [ 180/Pi*Atan2[Im[{Iz}],Re[{Iz}]] ]  ; In DomainZt_Cir ; }
        } 
	  }
	  

			
	  { Name Lac ;
      Value {
	  Term { Type Global; [ Im[{Uz}/{Iz}]/(2*Pi*Freq) ] ; In DomainZt_Cir ;  } } }		  

	EndIf  

	  
	  
	  
 
  EndIf  
	  
      // Not a good idea to use the functions...
  If (Irms_sec != 0)
      { Name Inductance_from_Flux ; Value { Term { Type Global; [ $Flux * 1e3/IA_sec0] ; In DomainDummy ; } } }
      { Name Inductance_from_MagEnergy ; Value { Term { Type Global; [ /*2 **/ $MagEnergy * 1e3/(Irms_sec*Irms_sec) ] ; In DomainDummy ; } } }
      { Name Resistance ; Value { Term { Type Global; [  $JouleLossesSec /(Irms_sec*Irms_sec) ] ; In DomainDummy ; } } }		
  EndIf

  If (Irms_pri != 0)
      { Name Inductance_from_Flux ; Value { Term { Type Global; [ $Flux * 1e3/IA_pri] ; In DomainDummy ; } } }
      { Name Inductance_from_MagEnergy ; Value { Term { Type Global; [ /*2 **/ $MagEnergy * 1e3/(Irms_pri*Irms_pri) ] ; In DomainDummy ; } } }
      { Name Resistance ; Value { Term { Type Global; [  $JouleLossesPri /(Irms_pri*Irms_pri) ] ; In DomainDummy ; } } }
  EndIf

  
  
    }
  }
}

 
//========================================================================================

PostOperation Get_LocalFields UsingPost MagStaDyn_av_js0_3D {
  Print[ a,  OnElementsOf Domain,  File StrCat[Dir, "a", ExtGmsh], LastTimeStepOnly ] ;
   If(_winding_model==STRANDED && _divJ_zero==DIVJ0_NONE)
     Print[ js0, OnElementsOf DomainS, File StrCat[Dir, "js0", ExtGmsh], LastTimeStepOnly ] ;
   EndIf
   If(_winding_model==STRANDED && _divJ_zero==DIVJ0_STRONG)
     Print[ hs, OnElementsOf DomainB, File StrCat[Dir, "hs", ExtGmsh], LastTimeStepOnly ] ;
     Print[ js, OnElementsOf DomainB, File StrCat[Dir, "js", ExtGmsh], LastTimeStepOnly ] ;
   EndIf
   If(_winding_model==STRANDED && _divJ_zero==DIVJ0_WEAK)
     Print[ xis,  OnElementsOf DomainS, File StrCat[Dir, "xis", ExtGmsh], LastTimeStepOnly ] ;
     Print[ dxis, OnElementsOf DomainS, File StrCat[Dir, "dxis", ExtGmsh], LastTimeStepOnly ] ;
     Print[ js0_dxis, OnElementsOf DomainS, File StrCat[Dir, "js0_corrected", ExtGmsh], LastTimeStepOnly ] ;
   EndIf
   If(_winding_model==MASSIVE)
     Print[ j, OnElementsOf DomainC,  File StrCat[Dir, "j", ExtGmsh], LastTimeStepOnly ] ;
   EndIf
  Print[ b,  OnElementsOf Domain,  File StrCat[Dir, "b",ExtGmsh], LastTimeStepOnly ] ;
 }


PostOperation Get_GlobalQuantities UsingPost MagStaDyn_av_js0_3D {
  Print[ Flux[DomainS], OnGlobal, Format TimeTable,
     File > StrCat[Dir,"Flux",ExtGnuplot], LastTimeStepOnly, StoreInVariable $Flux,
     SendToServer StrCat[po,"40Flux [Wb]"],  Color "LightYellow" ];

  Print[ Inductance_from_Flux, OnRegion DomainDummy, Format Table, LastTimeStepOnly,
     File StrCat[Dir,"InductanceF",ExtGnuplot],
     SendToServer StrCat[po,"50Inductance from Flux [mH]"], Color "LightYellow" ];

  Print[ MagEnergy[Domain], OnGlobal, Format TimeTable,
     File > StrCat[Dir,"ME",ExtGnuplot], LastTimeStepOnly, StoreInVariable $MagEnergy,
     SendToServer StrCat[po,"41Magnetic Energy [W]"],  Color "LightYellow" ];

  Print[ Inductance_from_MagEnergy, OnRegion DomainDummy, Format Table, LastTimeStepOnly,
     File StrCat[Dir,"InductanceE",ExtGnuplot],
     SendToServer StrCat[po,"51Inductance from Magnetic Energy [mH]"], Color "LightYellow" ];
  
  If (!_flag_circuit_coupling)
  
    Print[ JouleLosses[Primary], OnGlobal, Format TimeTable,
      File > StrCat[Dir,"WindingLoss",ExtGnuplot], LastTimeStepOnly, StoreInVariable $JouleLossesPri,
      SendToServer StrCat[po,"61Joule loss (pri) [W]"],  Color "LightYellow" ];
	 
    Print[ JouleLosses[Secondary], OnGlobal, Format TimeTable,
      File > StrCat[Dir,"WindingLoss",ExtGnuplot], LastTimeStepOnly, StoreInVariable $JouleLossesSec,
      SendToServer StrCat[po,"61Joule loss (sec) [W]"],  Color "LightYellow" ];
	
    Print[ Resistance, OnRegion DomainDummy, Format Table, LastTimeStepOnly,
      File StrCat[Dir,"ResistanceJ",ExtGnuplot],
      SendToServer StrCat[po,"71Winding resistance [Ohm]"], Color "LightYellow" ];
  
  EndIf
  
  If (_flag_circuit_coupling)
	If (IA_pri != 0)
      Print[ Rac, OnRegion VI_source_pri, Format TimeTable, File StrCat[Dir,"Rac_pri",ExtGnuplot],
        SendToServer StrCat[pcp,"4Rac_pri [Ohm]"]{0}, Color "LightYellow"];
	EndIf		

	If (IA_sec0 != 0 && IA_sec1 != 0)	
      Print[ Rac, OnRegion VI_source_sec, Format TimeTable, File StrCat[Dir,"Rac_sec",ExtGnuplot],
        SendToServer StrCat[pcs,"4Rac_sec [Ohm]"]{0}, Color "LightYellow"];
	EndIf	  
	  
	If (_analysis_type == STA)
      Print[ I, OnRegion VI_source_pri, Format TimeTable, File StrCat[Dir,"Current",ExtGnuplot],
        SendToServer StrCat[po,"I pri [A]"]{0}, Color "LightYellow"];
      Print[ U, OnRegion VI_source_pri, Format TimeTable, File StrCat[Dir,"Voltage",ExtGnuplot],
        SendToServer StrCat[po,"V pri [V]"]{0}, Color "LightYellow"];
      Print[ I, OnRegion Secondary0, Format TimeTable, File StrCat[Dir,"Current",ExtGnuplot],
        SendToServer StrCat[po,"I sec [A]"]{0}, Color "LightYellow"];
      Print[ U, OnRegion VI_source_sec, Format TimeTable, File StrCat[Dir,"Voltage",ExtGnuplot],
        SendToServer StrCat[po,"V sec [V]"]{0}, Color "LightYellow"];


		
	EndIf
	
    If (_analysis_type == DYN)
      Print[ Ipk, OnRegion VI_source_pri, Format TimeTable, File StrCat[Dir,"CurrentPriPk",ExtGnuplot],
        SendToServer StrCat[pcp,"0Ipk pri [A]"]{0}, Color "LightYellow"];
      Print[ Iph, OnRegion VI_source_pri, Format TimeTable, File StrCat[Dir,"CurrentPriPh",ExtGnuplot],
        SendToServer StrCat[pcp,"1Phase I pri [degree]"]{0}, Color "LightYellow"];
      Print[ Upk, OnRegion VI_source_pri, Format TimeTable, File StrCat[Dir,"VoltagePriPk",ExtGnuplot],
        SendToServer StrCat[pcp,"2Vpk pri [V]"]{0}, Color "LightYellow"];
      Print[ Uph, OnRegion VI_source_pri, Format TimeTable, File StrCat[Dir,"VoltagePriPh",ExtGnuplot],
        SendToServer StrCat[pcp,"3Phase V pri [degree]"]{0}, Color "LightYellow"];
		
      If (IA_pri != 0)
        Print[ Lac, OnRegion VI_source_pri, Format TimeTable, File StrCat[Dir,"Lac_pri",ExtGnuplot],
          SendToServer StrCat[pcp,"5Lac_pri [H]"]{0}, Color "LightYellow"];		
	  EndIf
	  
      Print[ Ipk, OnRegion VI_source_sec, Format TimeTable, File StrCat[Dir,"CurrentSecPk",ExtGnuplot],
        SendToServer StrCat[pcs,"0Ipk sec [A]"]{0}, Color "LightYellow"];
      Print[ Iph, OnRegion Secondary0, Format TimeTable, File StrCat[Dir,"CurrentSecPh",ExtGnuplot],
        SendToServer StrCat[pcs,"1Phase I sec [degree]"]{0}, Color "LightYellow"];
      Print[ Upk, OnRegion VI_source_sec, Format TimeTable, File StrCat[Dir,"VoltageSecPk",ExtGnuplot],
        SendToServer StrCat[pcs,"2Vpk sec [V]"]{0}, Color "LightYellow"];
      Print[ Uph, OnRegion VI_source_sec, Format TimeTable, File StrCat[Dir,"VoltageSecPh",ExtGnuplot],
        SendToServer StrCat[pcs,"3Phase V sec [degree]"]{0}, Color "LightYellow"];

	  If (IA_sec0 != 0 && IA_sec1 != 0)	
        Print[ Lac, OnRegion VI_source_sec, Format TimeTable, File StrCat[Dir,"Lac_sec",ExtGnuplot],
          SendToServer StrCat[pcs,"5Lac_sec [H]"]{0}, Color "LightYellow"];
	  EndIf
		
	EndIf	
		
  EndIf
}

DefineConstant[
  R_ = {"Analysis", Name "GetDP/1ResolutionChoices", Visible 1},
  C_ = {"-solve -v 3 -v2", Name "GetDP/9ComputeCommand", Visible 1},
  P_ = {"", Name "GetDP/2PostOperationChoices", Visible 1}
];