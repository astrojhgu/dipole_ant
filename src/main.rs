extern crate dipole_ant;

use num::{
    traits::{
        FloatConst
    }
};

use healpix_fits::{
    write_map
};

use dipole_ant::{
    free_dipole::{
        e_theta
    }
};

use scorus::{
    healpix::{
        pix::{
           pix2ang_ring 
        }
        , utils::{
            nside2npix
        }
        , rotation::{
            rotate_ring_pol
            , get_rotation_matrix
        }
    }
};

fn main() {
    let nside=128;
    let npix=nside2npix(nside);
    let mut stokes_I=Vec::new();
    let mut stokes_Q=Vec::new();
    let mut stokes_U=Vec::new();

    for (i,q,u) in (0..npix).map(|i| {
        let dir=pix2ang_ring(nside, i);
        let theta=dir.pol;
        let exx=e_theta(theta, 1.5).powi(2);
        let eyy=0.0;
        let stokes_I=exx+eyy;
        let stokes_Q=exx-eyy;
        let stokes_U=exx.sqrt()*eyy.sqrt()*1.0;
        (stokes_I, stokes_Q, stokes_U)
    }){
        stokes_I.push(i);
        stokes_Q.push(q);
        stokes_U.push(u);
    }
    write_map("i.fits", &[&stokes_I], false, true);
    write_map("q.fits", &[&stokes_Q], false, true);
    write_map("u.fits", &[&stokes_U], false, true);

    let mat=get_rotation_matrix(0.0, 45.0_f64.to_radians(), 0.0);
    let (stokes_I, stokes_Q, stokes_U)=rotate_ring_pol(&stokes_I, &stokes_Q, &stokes_U, &mat);

    write_map("i1.fits", &[&stokes_I], false, true);
    write_map("q1.fits", &[&stokes_Q], false, true);
    write_map("u1.fits", &[&stokes_U], false, true);

}
