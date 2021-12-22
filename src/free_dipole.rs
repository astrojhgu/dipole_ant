use num::{
    traits::{
        FloatConst
    }
    , Zero
};


type Ftype=f64;


pub fn e_theta(theta:Ftype, l_over_lambda: f64)->Ftype{
    let kl=2.0*Ftype::PI()*l_over_lambda;
    let st=theta.sin();
    if st==Ftype::zero(){
        Ftype::zero()
    }else{
        ((kl/2.0*theta.cos()).cos()-(kl/2.0).cos())/st
    }    
}
