ppw = "2Model parameters/3Wire parameters" ;


thick_insul00 = 5e-5;
Np = 7; //4
Ns = 6;

mm=1e-3;
ro  = (74-2.5)/2*mm; // outer-most radius
ri1 = (57.5+1.8)/2*mm + 1*mm;  // second-most radius //add to mesh easily
ri2 = (29.5-1)/2*mm; // center-post leg's radius
ri3 = (5.4+0.3)/2*mm; // inner-most hole

ro_clip = (4+1)/2*mm;


h  = (40.7+0.8)*mm; // bobbin height
ho = 59*mm; // Height

ag=3.8*mm;

rot_angle = Pi/2*0 - 0.1 + 1*Atan[ro*Tan[Pi/3]/(ro + ri2)];



rp00 = 0.0043/2;  
interwire_pri = 0.0035/4;//0.0035, 0.005065; // 0.003;// modified because the cross section shows 5 circles instead of 4 circles
rs00 = 0.00252/2; interwire_sec = 0.003; //0.003676; //0.003;//
strand_dia_pri_00 = 0.2*mm;
strand_dia_sec_00 = 0.2*mm;
strand_number_pri_00 = 200;
strand_number_sec_00 = 100;

z0 = 2*ro;
lb = z0; // length of the box
hb = 2*ho; // height of the box

// initial coordinate of the center of the coil in XY plane, where Y is height
xp0 = 0.022775; yp0 = 0.0148+0.001; zp0 = z0;
xs0 = 0.01781; ys0 = 0.01686; zs0 = z0;
xs1 = 0.02774; ys1 = 0.01686; zs1 = z0;


DefineConstant[
	rp = { rp00, Name StrCat[ppw, "/03Bundle radius of a primary bundle [m]" ], Highlight "LightGreen", ReadOnly 1 } 
	rs = { rs00, Name StrCat[ppw, "/04Bundle radius of a secondary bundle [m]" ], Highlight "LightGreen", ReadOnly 1 } 

	thick_insul = { thick_insul00, Name StrCat[ppw, "/09Insulation thickness [m]" ], Highlight "LightGreen", ReadOnly 1 }


	//Rpri = { (1/sigma_wire)*(Lz*2*/*(start_winding_area + 
];

// GetDP parameters

sigma_cu = 6e7;
SymmetryFactor = 1;

// Physical numbers...
// Volume
PRIMARY    = 1000; SKIN_PRIMARY    = 1e6 + PRIMARY;
SECONDARY0 = 2000; SKIN_SECONDARY0 = 1e6 + SECONDARY0;
SECONDARY1 = 2100; SKIN_SECONDARY1 = 1e6 + SECONDARY1;
CORE = 3000;

AIR = 5000;

// Surface
IN_PRI  = 11111;
OUT_PRI = 11112;
IN_SEC0  = 11113;
OUT_SEC0 = 11114;
IN_SEC1  = 11115;
OUT_SEC1 = 11116;

SURF_AIROUT = 11121;
