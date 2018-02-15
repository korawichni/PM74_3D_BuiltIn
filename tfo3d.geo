Include "tfo3d_data.pro";
Geometry.NumSubEdges = 100; // nicer display of curve

// Mesh control parameters
lc_wind = 1.5*rs; //2
lc_core = 2*lc_wind;
lc_inf  = 4*lc_core;
lc_ag = ag/4;

// Extrusion parameters
v_extrude = { -(h-ag)/2, -ag, -(h-ag)/2, -(ho-h)/2 };
nn_side = 6;
nn_cen = 4; // layer of the extrusion of the central middle part
nn_ag = 4; // layer of the airgap
nn_thick = 4; // layer of the base
nn_cen = { nn_thick, nn_ag, nn_cen, nn_thick };

// Layer of extruded wire
layer_pri = 40;
layer_sec0 = 40;
layer_sec1 = 40;

If (1)
	
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
aux_cen[] = Extrude {0,-(ho-h)/2,0} { Surface{surf_hole[],surf_cen[],surf_window[],surf_side[]}; Layers{nn_thick};  };
vol_airtube[] = aux_cen[{1:16:5}];
vol_core_cen[] = aux_cen[{21:39:6}];
vol_core_thick[] = aux_cen[{45:63:6}];
vol_core_side[] = aux_cen[{69:90:7}];
_surf_hole[] = aux_cen[{0:15:5}];
_surf_hole[] += aux_cen[{20:38:6}];
_surf_side[] = aux_cen[{68:89:7}];

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

// Extrude the side legs
aux_core_side[] = Extrude {0,-h,0} { Surface{_surf_side[]}; Layers{nn_side}; }; 
vol_core_side[] += aux_core_side[{1:22:7}];
_surf_side[] = aux_core_side[{0:21:7}];

// Extrude the base

aux_surf_base[] = Translate {0,-(ho+h)/2,0} { Duplicata { Surface{surf_window[]}; } };
aux_core_base = Extrude {0,-(ho-h)/2,0} { Surface{aux_surf_base[],_surf_side[]}; Layers{nn_thick}; };

vol_core_thick[] += aux_core_base[{1:19:6}];
vol_core_side[] += aux_core_base[{25:46:7}];


// Rotate the core
all_vol() = Volume '*';
Rotate {{0, 1, 0}, {0, 0, 0}, rot_angle} {Volume{all_vol()}; }

vol_core() = {vol_core_cen(),vol_core_side(),vol_core_thick()};


Printf('aux_core_side=',aux_core_base[]);


EndIf

//==================== Winding ============================//
If (1)

Macro DrawWinding

  vol_coil[] = {}; pnt0[] = {}; lncir[] = {};
  
  pnt0[] += newp; Point(newp) = {x0,y0,z0,lc0}; // center
  pnt0[] += newp; Point(newp) = {x0 + r,y0,z0,lc0};
  pnt0[] += newp; Point(newp) = {x0,y0  + r,z0,lc0};
  pnt0[] += newp; Point(newp) = {x0 - r,y0,z0,lc0};
  pnt0[] += newp; Point(newp) = {x0,y0 - r,z0,lc0};
  
  lncir[] += newl; Circle(newl) = {pnt0[1],pnt0[0],pnt0[2]};
  lncir[] += newl; Circle(newl) = {pnt0[2],pnt0[0],pnt0[3]};
  lncir[] += newl; Circle(newl) = {pnt0[3],pnt0[0],pnt0[4]};
  lncir[] += newl; Circle(newl) = {pnt0[4],pnt0[0],pnt0[1]};
  llcir = newll; Line Loop(newll) = lncir[];
  surf_cir = news; Plane Surface(news) = newll-1;
  
// Starting straight part: cable(0)
  aux() = Extrude{0.,0.,-z0}{ Surface{surf_cir}; Layers{_use_layers}; };
  vol_coil() += aux(1); // Get the volume from aux(1)
  surf_inout() = surf_cir; // Collect the surface
// Spiral part every one quarter for smooth curve
  For j In {1:4*Nturns-1}
	aux() = Extrude { {0,-(interwire+2*(r+thick_insul))/4,0}, {0,1,0} , {0,0,0} , 2*Pi/4 } {
	Surface{aux(0)}; Layers{_use_layers/2}; };  // Recombine;
	vol_coil() += aux(1);
  EndFor
  
// Ending straight part
  aux()   = Extrude{z0,0.,0}{ Surface{aux(0)}; Layers{_use_layers}; };
  vol_coil() += aux(1); // Get the volume from aux(1)
  surf_inout() += aux(0); // Collect the out surface
  
Return

//================ Primary winding ===========================//
x0 = xp0; y0 = yp0; z0 = lb; r = rp; lc0 = lc_wind;
Nturns = Np; interwire = interwire_pri;
_use_layers = layer_pri;
Call DrawWinding;
vol_coil_pri() = vol_coil();
surf_pri_inout() = surf_inout();
ll_pri_in = newll; Line Loop(newll) = Abs(Boundary{ Surface{surf_pri_inout(0)}; });
ll_pri_out = newll; Line Loop(newll) = Abs(Boundary{ Surface{surf_pri_inout(1)}; });
skin_pri[] = Abs(CombinedBoundary{ Volume{vol_coil_pri()}; });


skin_pri[] -= surf_pri_inout[];
Printf("surf_pri_inout() =", surf_pri_inout());
Printf("skin_pri() =", skin_pri());
//================ Secondary winding 0 ===========================//
x0 = xs0; y0 = ys0; z0 = lb; r = rs; lc0 = lc_wind;
Nturns = Ns; interwire = interwire_sec;
_use_layers = layer_sec0;
Call DrawWinding;
vol_coil_sec0() = vol_coil();
surf_sec0_inout() = surf_inout();
ll_sec0_in = newll; Line Loop(newll) = Abs(Boundary{ Surface{surf_sec0_inout(0)}; } );
ll_sec0_out = newll; Line Loop(newll) = Abs(Boundary{ Surface{surf_sec0_inout(1)}; } );

skin_sec0[] = Abs(CombinedBoundary{ Volume{vol_coil_sec0()}; });
skin_sec0[] -= surf_sec0_inout[];

//================ Secondary winding 1 ===========================//
x0 = xs1; y0 = ys1; z0 = lb; r = rs; lc0 = lc_wind;
Nturns = Ns; interwire = interwire_sec;
_use_layers = layer_sec1;
Call DrawWinding;
vol_coil_sec1() = vol_coil();
surf_sec1_inout() = surf_inout();
ll_sec1_in = newll; Line Loop(newll) = Abs(Boundary{ Surface{surf_sec1_inout(0)}; } );
ll_sec1_out = newll; Line Loop(newll) = Abs(Boundary{ Surface{surf_sec1_inout(1)}; } );

skin_sec1[] = Abs(CombinedBoundary{ Volume{vol_coil_sec1()}; });
skin_sec1[] -= surf_sec1_inout[];

EndIf

//================= Air Around ===================================//
If (1)
	
// Create a box
p = newp;
Point(p) = {lb,-hb/2,lb, lc_inf};
Point(p+1) = {lb,hb/2,lb, lc_inf};
Point(p+2) = {-lb,hb/2,lb, lc_inf};
Point(p+3) = {-lb,-hb/2,lb, lc_inf};

Point(p+4) = {lb,-hb/2,-lb, lc_inf};
Point(p+5) = {lb,hb/2,-lb, lc_inf};
Point(p+6) = {-lb,hb/2,-lb, lc_inf};
Point(p+7) = {-lb,-hb/2,-lb, lc_inf};

l = newl;
Line(l) = {p,p+1};
Line(l+1) = {p+1,p+2};
Line(l+2) = {p+2,p+3};
Line(l+3) = {p+3,p};
Line(l+4) = {p+4,p+5};
Line(l+5) = {p+5,p+6};
Line(l+6) = {p+6,p+7};
Line(l+7) = {p+7,p+4};
Line(l+8) = {p, p+4};
Line(l+9) = {p+1, p+5};
Line(l+10) = {p+2, p+6};
Line(l+11) = {p+3, p+7};

llb = newll;
Line Loop(llb) = {l+8, -(l+7), -(l+11), l+3}; //bottom
Line Loop(llb+1) = {l+9, l+5, -(l+10), -(l+1)}; //top
Line Loop(llb+2) = {l,l+1,l+2,l+3}; // face
Line Loop(llb+3) = {l+4,l+5,l+6,l+7}; // rear
Line Loop(llb+4) = {l+2, l+11, -(l+6), -(l+10)}; // left
Line Loop(llb+5) = {l, l+9, -(l+4), -(l+8)}; // right

surf_box[] += news; Plane Surface(news) = {llb};
surf_box[] += news; Plane Surface(news) = {llb+1};
surf_box[] += news; Plane Surface(news) = {llb+2,ll_pri_in,ll_sec0_in,ll_sec1_in};
surf_box[] += news; Plane Surface(news) = {llb+3};
surf_box[] += news; Plane Surface(news) = {llb+4};
surf_box[] += news; Plane Surface(news) = {llb+5,ll_pri_out,ll_sec0_out,ll_sec1_out};

sl_box = newsl; Surface Loop(sl_box) = surf_box[];

// Prepare surface loops to subtract from the box
sl_skin_pri = newsl; Surface Loop(sl_skin_pri) = skin_pri[];
sl_skin_sec0 = newsl; Surface Loop(sl_skin_sec0) = skin_sec0[];
sl_skin_sec1 = newsl; Surface Loop(sl_skin_sec1) = skin_sec1[];

vol_aux = {vol_core[],vol_ag[],vol_airtube[]};
skin_aux[] =  Abs(CombinedBoundary{ Volume{vol_aux[]}; });
sl_aux[] += newsl; Surface Loop(newsl) = Abs(CombinedBoundary{ Volume{vol_aux[]}; });

// Create air around
air_around = newv;
Volume(air_around) = {sl_box,sl_skin_pri,sl_skin_sec0,sl_skin_sec1,sl_aux};

EndIf

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

Physical Surface ("skin primary",1000) = skin_pri[];
Physical Surface ("skin secondary0",1001) = skin_sec0[];
Physical Surface ("skin aux",1002) = skin_aux[];
Physical Surface ("in pri",1003) = surf_pri_inout[0];
Physical Surface ("in sec0",1004) = surf_sec0_inout[0];
Physical Surface ("in sec1",1005) = surf_sec1_inout[0];
Physical Surface ("out pri",1006) = surf_pri_inout[1];
Physical Surface ("out sec0",1007) = surf_sec0_inout[1];
Physical Surface ("out sec1",1008) = surf_sec1_inout[1];

Physical Volume ("side legs",2000) = vol_core_side[];
Physical Volume ("center leg",3000) = vol_core_cen[];
Physical Volume ("air tube",4000) = vol_airtube[];
Physical Volume ("core thick",5000) = vol_core_thick[];
Physical Volume ("airgap",6000) = vol_ag[];
Physical Volume ("primary",7000) = vol_coil_pri[];
Physical Volume ("secondary0",8000) = vol_coil_sec0[];
Physical Volume ("secondary1",8100) = vol_coil_sec1[];
Physical Volume ("core",9000) = vol_core[];
Physical Volume ("air around",10000) = air_around;


// For aestetics
Recursive Color SkyBlue { Volume{air_around,vol_ag(),vol_airtube()};}
Recursive Color SteelBlue { Volume{vol_core()};}
Recursive Color Red  { Volume{vol_coil_pri()};}
Recursive Color Green{ Volume{vol_coil_sec0()};}
Recursive Color Yellow { Volume{vol_coil_sec1()};}