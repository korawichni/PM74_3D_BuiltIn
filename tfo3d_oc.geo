// geometry of winding and the mesh are ok

Include "tfo3d_data.geo";

Geometry.NumSubEdges = 100; // nicer display of curve

mdf= 2;

lc_wind = mdf*rs;
lc_core = 2*lc_wind;
lc_inf  = 4*lc_core;
lc_ag   = mdf*ag/4;


DefineConstant[
  _use_layers = {0, Choices{0,1}, Name "Geometry param./layers in extrusion?", ReadOnly 1}
  nls = {20, Name "Geometry param./number of layers"}
];


//=============================================================//

SetFactory("OpenCASCADE");


Macro DrawWinding
  // clear arrays
  llspiral() = {}; cable() = {}; cables() = {};

  // Need detail of the point {x0,y0,z0,lc0} and the radius
  // Need the NoOfTurn
  /*
  lc0  = r;
  // Straight part: cable(0)
  pnt0[] = {};
  pnt0[] += newp; Point(newp) = {x0,y0,z0, lc0}; // center
  pnt0[] += newp; Point(newp) = {x0 + r,y0,z0,lc0};
  pnt0[] += newp; Point(newp) = {x0,y0  + r,z0,lc0};
  pnt0[] += newp; Point(newp) = {x0 - r,y0,z0,lc0};
  pnt0[] += newp; Point(newp) = {x0,y0 - r,z0,lc0};

  lncir[] = {};
  lncir[] += newl; Circle(newl) = {pnt0[1],pnt0[0],pnt0[2]};
  lncir[] += newl; Circle(newl) = {pnt0[2],pnt0[0],pnt0[3]};
  lncir[] += newl; Circle(newl) = {pnt0[3],pnt0[0],pnt0[4]};
  lncir[] += newl; Circle(newl) = {pnt0[4],pnt0[0],pnt0[1]};
  llcir = newll; Line Loop(newll) = lncir[];
  surf_init(0) =news; Plane Surface(news) = llcir;
  */
  surf_init(0) =news; Disk(news) = {x0, y0, z0, r};

  aux()   = Extrude{0.,0.,-z0}{ Surface{surf_init(0)}; Layers{_use_layers*nls}; };
  cable() += aux(1); // Get the volume from aux(1)
  //Printf('aux=',aux());

  // Top extruded surface is aux(0)
  aux_l() = Boundary{ Surface{aux(0)}; }; // Get the line loop
  llspiral(0) = newll; Line Loop(newll) = aux_l(); // preparing to use thrusection

  For i In {1:section-1}
    aux_l() = Translate{0,-(interwire+r*2+thick_insul*2)*NoOfTurn/(section-1),0}{ Duplicata{ Line{aux_l()}; } };
    Rotate {{0, 1, 0}, {0, 0, 0}, angle/(section-1)} { Line{aux_l()}; }
    llspiral() += newll; Line Loop(newll) = aux_l(); // collect the line loop to use thrusection
  EndFor

  // Spiral part
  For i In {1:#llspiral()-nn+1:nn}
    k =  (i < #llspiral()-nn) ? 0 : 1 ;
    cable() += newv; ThruSections(newv) = llspiral({i-1:i-1+nn-k});
  EndFor

  /*
  For i In {1:#llspiral()-2+1:2}
    k =  (i < #llspiral()-2) ? 0 : 1 ;
    cable() += newv; ThruSections(newv) = llspiral({i-1:i-1+2-k});
  EndFor
  */
  // Ending straight part
  surf_init(1) = news; Plane Surface(news) = llspiral(#llspiral()-1);
  aux2() = Extrude{z0,0.,0.}{ Surface{surf_init(1)}; Layers{_use_layers*nls}; };
  cable() += aux2(1);

  vcable() = BooleanFragments{ Volume{cable()}; Delete; }{}; // be careful of cables overlapping with each other
Return

//======================================================================================

// The number of section cannot be too low otherwise the spiral will deviate from the path and hit the core

nnp = 10; // for some reason, it cannot be 10
nns = 12;
//================ Primary winding ===========================//
x0 = xp0; y0 = yp0; z0 = zp0; r = rp;
section = nnp*Np; NoOfTurn = Np; interwire = interwire_pri; angle = Pi*(2*Np-0.5);
nn = nnp;
Call DrawWinding;
winding_pri() = vcable();


//================ Secondary winding 0 ===========================//
x0 = xs0; y0 = ys0; z0 = zs0; r = rs;
section = nns*Ns; NoOfTurn = Ns; interwire = interwire_sec; angle = Pi*(2*Ns-0.5);
nn = nns;

Call DrawWinding;
winding_sec0() = vcable();
//BooleanUnion{ Volume{winding_sec0()}; }{}

//================ Secondary winding 1 ===========================//
x0 = xs1; y0 = ys1; z0 = zs1; r = rs;
section = nns*Ns; NoOfTurn = Ns; interwire = interwire_sec; angle = Pi*(2*Ns-0.5);
nn = nns;

Call DrawWinding;
winding_sec1() = vcable();

Characteristic Length { PointsOf{ Volume{winding_sec0(),winding_sec1(),winding_pri()}; } } = lc_wind;


//========================== Core ==================================//

vC_out() +=newv; Cylinder(newv) = {0, -ho/2, 0, 0, ho, 0, ro, 2*Pi}; // Outer-most cylinder
vC_in=newv; Cylinder(newv) = {0, -h/2, 0, 0, h, 0, ri1, 2*Pi};  // Second-most cylinder

vC_out() = BooleanDifference{ Volume{vC_out()}; Delete;}{ Volume{vC_in}; Delete;};
vC_cen = newv; Cylinder(newv) = {0, -h/2, 0, 0, h, 0, ri2, 2*Pi}; // Center-post leg

// Auxiliary cylinders to cut off the core
aux_[0] = newv; Wedge(aux_[0]) = {0, 0, 0, ro*Tan(Pi/3), ro, ho, 0};
Rotate {{-1, 0, 0}, {0, 0, 0}, Pi/2} {Volume{aux_[0]}; }
Translate {0, -ho/2, ro + ri2} {Volume{aux_[0]}; }

aux_[] += Symmetry {1, 0, 0, 0} { Duplicata{ Volume{aux_[0]}; } };

rot_angle = Pi/2*0 - 0.1 + 1*Atan(ro*Tan(Pi/3)/(ro + ri2));
aux_[] += Symmetry {0, 0, 1, 0} { Duplicata{ Volume{aux_[]}; } };
Rotate {{0, 1, 0}, {0, 0, 0}, rot_angle} {Volume{aux_()}; }

// Hole for clipping
aux_1[0] = newv;
Cylinder(aux_1[0]) = {ro, ho/2, -0, 0, -ho, 0, ro_clip, 2*Pi};
aux_1[] += newv;
Cylinder(aux_1[1]) = {-ro, ho/2, -0, 0, -ho, 0, ro_clip, 2*Pi};
Rotate {{0, 1, 0}, {0, 0, 0}, rot_angle} {Volume{aux_1()}; }

// Cut off the core
vC_out() = BooleanDifference{ Volume{vC_out()}; Delete; }{ Volume{aux_(),aux_1[]}; Delete; }; // Volume{aux_()};

// Airgap cylinder
vAg=newv; Cylinder(newv) = {0, -ag/2, 0, 0, ag, 0, ri2, 2*Pi};
// Make the airgap
vC() = BooleanFragments{ Volume{vC_cen}; Delete; }{ Volume{vAg}; Delete; };


// Stack the outer part and the central leg
vC() = BooleanFragments{ Volume{vC_out()}; Delete; }{ Volume{vC()}; Delete; };

// Make the hole
//hole_cen =newv; Cylinder(newv) = {0, -ho/2, 0, 0, -ho, 0, ri3, 2*Pi}; // hole in the center
hole_cen_up =newv; Cylinder(newv) = {0, ho/2, 0, 0, -(ho-ag)/2, 0, ri3, 2*Pi}; // hole in the center
hole_cen_dn =newv; Cylinder(newv) = {0, -ho/2, 0, 0, (ho-ag)/2, 0, ri3, 2*Pi}; // hole in the center
// Cut off the hole
//vC() = BooleanDifference{ Volume{vC()}; Delete; }{Volume{hole_cen_up,hole_cen_dn}; Delete; };
vC() = BooleanFragments{ Volume{vC()}; Delete; }{Volume{hole_cen_up,hole_cen_dn}; Delete; };

// Printf("vC_cen", vC_cen()); //0:20 1:21 2:22 3:23

// Without this the mesh will be too big -> error of self intersecting surface
Characteristic Length { PointsOf{ Volume{vC()}; } } = lc_core;


//================= Air Around ===================================//
air_around = newv; Box(air_around) = {-zs0, -zs0, zs0, 2*zs0, 2*zs0, -2*zs0};
Characteristic Length { PointsOf{ Volume{air_around}; } } = lc_inf;

vol_model() = BooleanFragments{ Volume{winding_pri(),winding_sec0(),winding_sec1(),vC(),air_around}; Delete; }{};
vol_model() -= {winding_pri(),winding_sec0(),winding_sec1()};

nn = #vol_model()-1;
vol_core() = vol_model({1,4,6});
vol_model() -= vol_core();
vol_air() = vol_model();
Characteristic Length { PointsOf{ Volume{vol_air(0)}; } } = lc_ag; // Airgap

Printf("", vol_air());

//------------------
// Physical regions
//------------------
Physical Volume ("Primary", PRIMARY+0) = winding_pri();
Physical Volume ("Secondary 0",SECONDARY0+0) = winding_sec0();
Physical Volume ("Secondary 1",SECONDARY1+0) = winding_sec1();

/*
nnp = #winding_pri()-1;
Physical Volume ("Primary helix", PRIMARY+0) = winding_pri({1:nnp-1});
Physical Volume ("Primary in",    PRIMARY+1) = winding_pri(0);
Physical Volume ("Primary out",   PRIMARY+2) = winding_pri(nnp);

nns0 = #winding_sec0()-1;
Physical Volume ("Secondary 0 helix",SECONDARY0+0) = winding_sec0({1:nns0-1});
Physical Volume ("Secondary 0 in",   SECONDARY0+1) = winding_sec0(0);
Physical Volume ("Secondary 0 out",  SECONDARY0+2) = winding_sec0(nns0);

nns1 = #winding_sec1()-1;
Physical Volume ("Secondary 1 helix",SECONDARY1+0) = winding_sec1({1:nns1-1});
Physical Volume ("Secondary 1 in",   SECONDARY1+1) = winding_sec1(0);
Physical Volume ("Secondary 1 out",  SECONDARY1+2) = winding_sec1(nns1);
*/

Physical Volume ("Air", AIR) = vol_air();
Physical Volume ("Core", CORE) = vol_core();


// Manually fill in
Physical Surface ("In surface primary", IN_PRI) = {97};
Physical Surface ("Out surface primary", OUT_PRI) = {98};

Physical Surface ("In surface secondary0", IN_SEC0) = {99};
Physical Surface ("Out surface secondary0", OUT_SEC0) = {100};

Physical Surface ("In surface secondary1", IN_SEC1) = {101};
Physical Surface ("Outsurface secondary1", OUT_SEC1) = {102};

Physical Surface ("Air boundary", SURF_AIROUT) = {85,86,87,88,89,90};


bnd_pri()  = CombinedBoundary{Volume{winding_pri()};};  bnd_pri() -= {97,98};
bnd_sec0() = CombinedBoundary{Volume{winding_sec0()};}; bnd_sec0()-= {99,100};
bnd_sec1() = CombinedBoundary{Volume{winding_sec1()};}; bnd_sec1()-= {100,101};

Physical Surface ("skin primary", SKIN_PRIMARY) = bnd_pri();
Physical Surface ("skin secondary 0", SKIN_SECONDARY0) = bnd_sec0();
Physical Surface ("skin secondary 1", SKIN_SECONDARY1) = bnd_sec1();


// For aestetics
Recursive Color SkyBlue { Volume{vol_air()};}
Recursive Color SteelBlue { Volume{vol_core()};}

Recursive Color Red  { Volume{winding_pri()};}
Recursive Color Green{ Volume{winding_sec0()};}
Recursive Color Cyan { Volume{winding_sec1()};}


Cohomology(1) {{AIR, CORE},{}}; // Hphi formulation
