use num::{
    traits::{
        FloatConst
    }
};

use healpix_fits::write_map;

use scorus::{
    healpix::{
        interp::{
            get_interpol_ring
        }
        , utils::{
            nside2npix
        }
    }
    , coordinates::{
        SphCoord
    }
};

pub fn main(){
    let mut context=necrs::RawNecContext::new();
    let nsegs=vec![39,37,35,33,31,29,27,25,23,21,19,17,15,13,11];
    let coords=vec![
        vec![0.000000,-1.529613,0.000000,0.000000,1.529613,0.000000],
        vec![0.331851,-1.369949,0.000000,0.331851,1.369949,0.000000],
        vec![0.629056,-1.226972,0.000000,0.629056,1.226972,0.000000],
        vec![0.895248,-1.098906,0.000000,0.895248,1.098906,0.000000],
        vec![1.133653,-0.984199,0.000000,1.133653,0.984199,0.000000],
        vec![1.347191,-0.881482,0.000000,1.347191,0.881482,0.000000],
        vec![1.538427,-0.789457,0.000000,1.538427,0.789457,0.000000],
        vec![1.709699,-0.707060,0.000000,1.709699,0.707060,0.000000],
        vec![1.863090,-0.633273,0.000000,1.863090,0.633273,0.000000],
        vec![2.000479,-0.567157,0.000000,2.000479,0.567157,0.000000],
        vec![2.123516,-0.507975,0.000000,2.123516,0.507975,0.000000],
        vec![2.233727,-0.454939,0.000000,2.233727,0.454939,0.000000],
        vec![2.332431,-0.407467,0.000000,2.332431,0.407467,0.000000],
        vec![2.420823,-0.364922,0.000000,2.420823,0.364922,0.000000],
        vec![2.499995,-0.326847,0.000000,2.499995,0.326847,0.000000],
    ];

    let tl_connection=vec![
        vec![1,20,2,19],
        vec![2,19,3,18],
        vec![3,18,4,17],
        vec![4,17,5,16],
        vec![5,16,6,15],
        vec![6,15,7,14],
        vec![7,14,8,13],
        vec![8,13,9,12],
        vec![9,12,10,11],
        vec![10,11,11,10],
        vec![11,10,12,9],
        vec![12,9,13,8],
        vec![13,8,14,7],
        vec![14,7,15,6],
    ];
    let dia=5e-3;
    let tl_impedance=-73.0;
    for (i, (coord, &n)) in coords.iter().zip(nsegs.iter()).enumerate(){
        context.nec_wire(i as i32+1, n, coord[0],coord[1], coord[2], coord[3], coord[4], coord[5], dia/2.0, 1.0, 1.0);
    }
    context.nec_gm_card(0, 0, 0.0, -90.0, 90.0, 0.0,  0.0, 0.5, 0);
    context.nec_geometry_complete(1);
    context.nec_gn_card(1, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    context.nec_ex_card(0, 15, 6, 0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    for ids in &tl_connection{
        let t1=ids[0];
        let s1=ids[1];
        let t2=ids[2];
        let s2=ids[3];
        context.nec_tl_card(t1, s1, t2, s2, tl_impedance, 0.0, 0.0, 0.0, 0.0, 0.0);
    }

    context.nec_fr_card(0, 1, 60.0, 0.0);
    


    let nside=64;
    let npix=nside2npix(nside);
    let angular_resolution=(4.0*f64::PI()/npix as f64).sqrt().to_degrees();
    println!("{}", angular_resolution);
    let nphi=(360.0/angular_resolution*4.0) as usize+1;
    let ntheta=(nphi-1)/2+1;
    let dphi=360.0/(nphi-1) as f64;
    let dtheta=180.0/(ntheta-1) as f64;
    
    //context.nec_rp_card(0, ntheta as i32, nphi as i32, 1, 0, 0, 0, 0.0, 0.0, dtheta, dphi, 0.0, 0.0);
    let (thetas, phis)=context.rp_from_npix(npix*16, 0, 1, 0, 0, 0, 0.0, 0.0);

    let mut data=vec![0.0; npix];
    let mut wgt=vec![0.0; npix];

    for (i, &theta) in thetas.iter().enumerate(){
        for (j, &phi) in phis.iter().enumerate(){
            let g=(context.nec_gain(0, i as i32, j as i32)/10.0).exp();
            let dir=SphCoord::new(theta.to_radians(), phi.to_radians());
            let (pix, w)=get_interpol_ring(nside, dir);
            for (&p,&w) in pix.iter().zip(w.iter()){
                wgt[p]+=w;
                data[p]+=w*g;
            }
        }
    }
    for (d,&w) in data.iter_mut().zip(wgt.iter()){
        if w>0.0{
            *d/=w;
        }
    }

    write_map("a.fits", &[&data], false, true);
}
