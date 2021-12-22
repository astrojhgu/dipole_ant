use clap::{
    App
    , Arg
};

use std::{
    fs::File
    , io::Write
};

use num::{
    traits::FloatConst
};

use scorus::{
    healpix::{
        nside2npix
        , pix2ang_ring
        , interp::natural_interp_ring
        , npix2nside
    }
};

use healpix_fits::{
    read_map
    , standard_column_names
};

fn main(){
    let matches=App::new("show pol")
    .arg(
        Arg::new("t.fits")
        .short('t')
        .long("t")
        .takes_value(true)
        .value_name("t.fits")
        .required(false)
    )
    .arg(
        Arg::new("q.fits")
        .short('q')
        .long("q")
        .takes_value(true)
        .value_name("q.fits")
        .required(false)
    )
    .arg(
        Arg::new("u.fits")
        .short('u')
        .long("u")
        .takes_value(true)
        .value_name("u.fits")
        .required(false)
    )
    .arg(
        Arg::new("pol_map")
        .short('i')
        .long("input")
        .takes_value(true)
        .value_name("pol map")
        .required(false)
    )
    .arg(
        Arg::new("outname")
        .short('o')
        .long("out")
        .takes_value(true)
        .value_name("out file name")
        .required(true)
    )
    .arg(
        Arg::new("sampling resolution")
        .short('s')
        .long("samp")
        .takes_value(true)
        .value_name("sampling resolution")
        .required(true)
    )
    .get_matches();


    let pol_map=matches.value_of("pol_map").map_or_else(||{
        let column_names=standard_column_names(1);
        let ref_column_names:Vec<_>=column_names.iter().map(|x| x.as_str()).collect();
        ["t.fits", "q.fits", "u.fits"].into_iter().map(|pname| {
            let fname=matches.value_of(pname).unwrap();
            let mut data=Vec::new();
            let mut result=read_map::<f64>(fname, &ref_column_names, 1);
            std::mem::swap(&mut data, &mut result[0]);
            data
        }).collect()
    }, |v|{
        let column_names=standard_column_names(3);
        let ref_column_names:Vec<_>=column_names.iter().map(|x| x.as_str()).collect();
        read_map::<f64>(v, &ref_column_names, 1)
    });



    let npix=pol_map[0].len();
    let nside=npix2nside(npix);


    let nside_samp=matches.value_of("sampling resolution").unwrap().parse::<usize>().unwrap();
    let npix_samp=nside2npix(nside_samp);

    //let t_max=pol_map[0].iter().zip(pol_map[1].iter()).map(|(&q, &u)| (q.powi(2)+u.powi(2)).sqrt()).fold(0.0, |a,b|{if a>b{a}else{b}});
    let t_max=pol_map[0].iter().fold(0.0, |a,&b|{if a>b{a}else{b}});
    let bar_unit_arcsec=(4.0*f64::PI()/npix_samp as f64).sqrt().to_degrees()*3600.0/2.0;

    let mut outfile=File::create(matches.value_of("outname").unwrap()).unwrap();
    writeln!(&mut outfile, "# Region file format: DS9 version 4.1").unwrap();
    writeln!(&mut outfile, "global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1").unwrap();
    writeln!(&mut outfile, "fk5").unwrap();
    for ipix in 0..npix_samp{
        let sph=pix2ang_ring::<f64>(nside_samp, ipix);
        let t=natural_interp_ring(nside, &pol_map[0], sph);
        let q=natural_interp_ring(nside, &pol_map[1], sph);
        let u=natural_interp_ring(nside, &pol_map[2], sph);
        let p=(q.powi(2)+u.powi(2)).sqrt();
        let psi=f64::atan2(u, q)/2.0;
        let angle=(f64::PI()/2.0-psi).to_degrees();
        let lat=(f64::PI()/2.0-sph.pol).to_degrees();
        let lon=sph.az.to_degrees();
        let bar_len=p/t_max*bar_unit_arcsec;
        let r=t/t_max*bar_unit_arcsec;
        writeln!(&mut outfile, "circle({},{},{}\"", lon, lat, r).unwrap();
        writeln!(&mut outfile, "ellipse({},{},{}\",1\",{})", lon, lat, bar_len, angle).unwrap();
    }
}
