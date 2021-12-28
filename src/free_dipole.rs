use num::{
    traits::{
        FloatConst
    }
    , Zero
};


type Ftype=f64;


pub fn e_theta(theta:Ftype, l_over_lambda: Ftype)->Ftype{
    let kl=2.0*Ftype::PI()*l_over_lambda;
    let st=theta.sin();
    if st==Ftype::zero(){
        Ftype::zero()
    }else{
        ((kl/2.0*theta.cos()).cos()-(kl/2.0).cos())/st
    }    
}

pub fn stokes(theta: Ftype, l_over_lambda: Ftype)->(Ftype, Ftype, Ftype){
    let exx=e_theta(theta, l_over_lambda).powi(2);
    let eyy=0.0;
    let stokes_I=exx+eyy;
    let stokes_Q=exx-eyy;
    let stokes_U=exx.sqrt()*eyy.sqrt()*1.0;
    (stokes_I, stokes_Q, stokes_U)
}