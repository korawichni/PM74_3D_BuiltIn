Include "tfo3d_data.pro";

lc_wind = 1.5*rs; //2
lc_core = 2*lc_wind;
lc_inf  = 4*lc_core;
lc_ag = ag/4;

nn_side = 6;
nn_cen = 4; // layer of the extrusion of the central leg
nn_ag = 4;
nn_thick = 4;

v_extrude = {-(h-ag)/2,-ag,-(h-ag)/2,-(ho-h)/2};
nn_cen = {nn_thick,nn_ag,nn_cen,nn_thick};

origin = newp;  Point(origin) = {0,    ho/2, 0,   lc_core};
pntH[] += newp; Point(newp)   = {-ri3, ho/2, 0,   lc_core};
pntH[] += newp; Point(newp)   = {0,    ho/2, ri3, lc_core};

pntC[] += newp; Point(newp)   = {-ri2, ho/2, 0,   lc_core};
pntC[] += newp; Point(newp)   = {0,    ho/2, ri2, lc_core};

pntR[] += newp; Point(newp)   = {-ri1, ho/2, 0,   lc_core};
pntR[] += newp; Point(newp)   = {-0.0181266,    ho/2, 0.0247154, lc_core};

pntO[] += newp; Point(newp)   = {-(ro-ro_clip), ho/2, 0,   lc_core};
pntO[] += newp; Point(newp)   = {-ro,    ho/2, 0, lc_core};
pntO[] += newp; Point(newp)   = {-ro,    ho/2, ro_clip, lc_core};
pntO[] += newp; Point(newp)   = {-0.0228868,    ho/2, 0.0274637, lc_core};



lnH[] += newl; Line(newl) = {origin,pntH[0]};
lnH[] += newl; Circle(newl) = {pntH[0],origin,pntH[1]};
lnH[] += newl; Line(newl) = {pntH[1],origin};

lnC[] += newl; Line(newl) = {pntH[0],pntC[0]};
lnC[] += newl; Circle(newl) = {pntC[0],origin,pntC[1]};
lnC[] += newl; Line(newl) = {pntC[1],pntH[1]};

lnR[] += newl; Line(newl) = {pntC[0],pntR[0]};
lnR[] += newl; Circle(newl) = {pntR[0],origin,pntR[1]};
lnR[] += newl; Line(newl) = {pntR[1],pntC[1]};

lnO[] += newl; Line(newl) = {pntR[0],pntO[0]};
lnO[] += newl; Circle(newl) = {pntO[0],pntO[1],pntO[2]};
lnO[] += newl; Circle(newl) = {pntO[2],origin,pntO[3]};
lnO[] += newl; Line(newl) = {pntO[3],pntR[1]};
// Center hole
llH[] += newll; Line Loop(newll) = lnH[];
surf_hole[] += news; Plane Surface(news) = llH[];
// Central leg
llC[] += newll; Line Loop(newll) = {lnC[],-lnH[1]};
surf_cen[] += news; Plane Surface(news) = llC[];
// Winding window
llR[] += newll; Line Loop(newll) = {lnR[],-lnC[1]};
surf_window[] += news; Plane Surface(news) = llR[];
// Side leg
llO[] += newll; Line Loop(newll) = {lnO[],-lnR[1]};
surf_side[] += news; Plane Surface(news) = llO[];

//================ Draw the top surface for extrusion ======================//
surf_aux[] = Symmetry {0,0,1,0} { Duplicata {Surface{surf_hole[],surf_cen[],surf_window[],surf_side[]}; } };
//surf_aux[] = -surf_aux[];
Printf("surf_aux=",surf_aux[]);
surf_hole[] += surf_aux[0];
surf_cen[] += surf_aux[1];
surf_window[] += surf_aux[2];
surf_side[] += surf_aux[3];

surf_aux[] = Symmetry {1,0,0,0} { Duplicata {Surface{surf_hole[],surf_cen[],surf_window[],surf_side[]}; } };
Printf("surf_aux=",surf_aux[]);
surf_hole[] += surf_aux[{0,1}];
surf_cen[] += surf_aux[{2,3}];
surf_window[] += surf_aux[{4,5}];
surf_side[] += surf_aux[{6,7}];

//================ Extrude the center leg ======================//
// Extrude the first thickness
aux_cen[] = Extrude {0,-(ho-h)/2,0} { Surface{surf_hole[],surf_cen[],surf_window[]}; Layers{nn_thick};  };
vol_airtube[] = aux_cen[{1:16:5}];
vol_core_cen[] = aux_cen[{21:39:6}];
vol_core_thick[] = aux_cen[{45:63:6}];
_surf_hole[] = aux_cen[{0:15:5}];
_surf_hole[] += aux_cen[{20:38:6}];
Printf("vol_aux_cen=",aux_cen[]);
Printf("_surf_hole=",_surf_hole[]);

//Physical Surface ("tmp",10001) = _surf_hole[];

// Extrude the center leg 4 steps 
For i In {0:3}
  aux_cen[] = Extrude {0,v_extrude[i],0} { Surface{_surf_hole[]}; Layers{nn_cen[i]};  };
  _surf_hole[] = aux_cen[{0:15:5}];
  _surf_hole[] += aux_cen[{20:38:6}];
  
  // Collect airgap
  If (i==1)
    vol_ag[] = aux_cen[{1:16:5}];
	vol_ag[] += aux_cen[{21:39:6}];
  Else
    vol_airtube[] += aux_cen[{1:16:5}];
    vol_core_cen[] += aux_cen[{21:39:6}];
  EndIf
  
EndFor


//==================== Winding ============================//
/* aux_cen[] = Extrude {0,-(h-ag)/2,0} { Surface{_surf_hole[]}; Layers{nn_cen[1]};  };
Printf("aux_cen2=",aux_cen[]);
_surf_hole[] = aux_cen[{0:15:5}];
_surf_hole[] += aux_cen[{20:38:6}];
vol_airtube[] += aux_cen[{1:16:5}];
vol_core_cen[] += aux_cen[{21:39:6}];

aux_cen[] = Extrude {0,-ag,0} { Surface{_surf_hole[]}; Layers{nn_airgap};  };
_surf_hole[] = aux_cen[{0:15:5}];
_surf_hole[] += aux_cen[{20:38:6}];
vol_airtube[] += aux_cen[{1:16:5}];
vol_core_cen[] += aux_cen[{21:39:6}];

aux_cen[] = Extrude {0,-(h-ag)/2,0} { Surface{_surf_hole[]}; Layers{nn_cen[1]};  };
_surf_hole[] = aux_cen[{0:15:5}];
_surf_hole[] += aux_cen[{20:38:6}];
vol_airtube[] += aux_cen[{1:16:5}];
vol_core_cen[] += aux_cen[{21:39:6}];

aux_cen[] = Extrude {0,-(ho-h)/2,0} { Surface{_surf_hole[]}; Layers{nn_cen[0]};  };
vol_airtube[] += aux_cen[{1:16:5}];
vol_core_cen[] += aux_cen[{21:39:6}];

*/

vol_core_thick[] += Translate {0,-(ho+h)/2,0} { Duplicata { Volume{vol_core_thick[]}; } };


vol_core_side[] = Extrude {0,-ho,0} { Surface{surf_side[]}; Layers{nn_side}; }; 



Physical Surface ("surface hole",1000) = surf_hole[];
Physical Surface ("surface center",1001) = surf_cen[];
Physical Surface ("surface window",1002) = surf_window[];
Physical Surface ("surface side",1003) = surf_side[];

Physical Volume ("side legs",2000) = vol_core_side[];
Physical Volume ("center leg",3000) = vol_core_cen[];
Physical Volume ("air tube",4000) = vol_airtube[];
Physical Volume ("core thick",5000) = vol_core_thick[];
Physical Volume ("airgap",6000) = vol_ag[];

// pntT[] += newp; Point(newp) = {0,ho/2,0,lc_core};
// pntT[] += newp; Point(newp) = {0,ho/2,0,lc_core};
// pntT[] += newp; Point(newp) = {0,ho/2,0,lc_core};
// pntT[] += newp; Point(newp) = {0,ho/2,0,lc_core};
// pntT[] += newp; Point(newp) = {0,ho/2,0,lc_core};