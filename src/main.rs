use std::f64::consts::PI;

fn wigner(beta: f64, omega: f64) -> f64 {
    let u = omega * beta;

    return 1.0 + u * u / 24.0;
}

fn bell(beta: f64, omega: f64) -> f64 {
    let u = omega * beta;

    return (0.5 * u / f64::sin(0.5 * u)).abs();
}

fn skodje_truhlar(beta: f64, omega: f64, v0: f64) -> f64 {
    let alpha = 2.0 * PI / omega;

    if beta > alpha {
        // When omega is large and low T
        let kappa = (((beta - alpha) * v0).exp() - 1.0) * beta / (beta - alpha);

        //println!("Skodje, beta*omega: {}", beta * omega);
        //println!("beta*omega > twopi");

        return kappa;
    } else if beta < alpha {
        // When omega is small and high T
        let dum0 = PI * beta / alpha;
        let dum1 = dum0 / dum0.sin();
        let dum2 = ((beta - alpha) * v0).exp() * beta / (beta - alpha);
        let kappa = dum1 + dum2;

        //println!("Skodje, beta*omega: {}", beta * omega);
        //println!("beta*omega < twopi");

        return kappa;
    } else {
        panic!("Wrong mode-argument in Skodje_Truhlar routine");
    }
}

fn skodje_truhlar_exact(beta: f64, omega: f64, v0: f64) -> f64 {
    let alpha = 2.0 * PI / omega;

    const NMAX: usize = 100;

    let mut res = 0.0;

    for n in 0..=NMAX {
        let numerator = 1.0 - ((beta - (n + 1) as f64 * alpha).exp() * v0);
        let denom = (n + 1) as f64 * alpha - beta;
        let dum = 1.0 / (n as f64 * alpha + beta);
        let mut alter = 1.0;

        if n % 2 == 1 {
            alter = -1.0;
        }

        res += alter * beta * (numerator / denom + dum);
    }

    return res;
}

fn eckart(beta: f64, omega: f64, vf: f64, vb: f64, de: f64, emax: f64) -> (f64, f64) {
    let n_emax = (emax / de) as usize;
    let mut to_int1 = vec![0.0; n_emax];
    let mut to_int2 = vec![0.0; n_emax];

    let alpha1 = 2.0 * PI * vf / omega;
    let alpha2 = 2.0 * PI * vb / omega;

    let a = (vf.sqrt() + vb.sqrt()).powi(2);
    let b = vf - vb;
    let d = ( (b*b - a*a).powi(2) / (a*a*a) / 8.0).sqrt() /  omega;

    let e0 = 0.0;
    for i in 0..n_emax {
        let e = e0 + (i as f64) * de;
        let csi = e / vf;

        let ptun1 = tunprop1(alpha1, alpha2, csi);
        let ptun2 = tunprop2(a, b, d, e);

        to_int1[i] = (-e * beta).exp() * ptun1;
        to_int2[i] = (-e * beta).exp() * ptun2;

    }
    let n_emin = 0;

    //let mut res1 = 0.0;
    //for i in n_emin+1..n_emax {
    //    res1 += 0.5*(to_int1[i]+to_int1[i-1])*de;
    //}

    let res1 = simpson_integrate(&to_int1, n_emin, n_emax, de);
    let kappa1 = res1 * (beta * vf).exp() * beta;

    let res2 = simpson_integrate(&to_int2, n_emin, n_emax, de);
    let kappa2 = res2 * (beta * vf).exp() * beta;

    return (kappa1, kappa2);
}

fn tunprop1(alpha1: f64, alpha2: f64, csi: f64) -> f64 {
    let denom = 1.0 / alpha1.sqrt() + 1.0 / alpha2.sqrt();
    let twopi_a = 2.0 * (alpha1 * csi).sqrt() / denom;
    let twopi_b = 2.0 * ((alpha1 * csi) + (alpha1 - alpha2).abs()).sqrt() / denom;
    let twopi_d = 2.0 * ((alpha1 * alpha2) - (PI * PI / 4.0)).sqrt();


    let csh_amb = (twopi_a - twopi_b).cosh();
    let csh_apb = (twopi_a + twopi_b).cosh();
    let csh_d = twopi_d.cosh();

    let res = 1.0 - (csh_amb + csh_d) / (csh_apb + csh_d);
    return res;
}

fn tunprop2(a: f64, b: f64, d: f64, e: f64) -> f64 {
    let twopi_d = 2.0 * PI * d;

    let alpha = twopi_d * (2.0 * e).sqrt();
    let beta  = twopi_d * (2.0 * (e - b)).sqrt();
    let delta = twopi_d * ((2.0 * a) - (1.0 / (2.0 * d * d))).sqrt();

    let csh_amb = (alpha - beta).cosh();
    let csh_apb = (alpha + beta).cosh();
    let csh_d = delta.cosh();

    let res = 1.0 - (csh_amb + csh_d) / (csh_apb + csh_d);
    return res;
}

fn simpson_integrate(func: &[f64], nmin: usize, nmax: usize, step: f64) -> f64 {
    let mut s0 = 0.0;
    let mut s1 = 0.0;
    let mut s2 = 0.0;

    if nmin > nmax {
        panic!("Wrong boundaries in Simpson Integration: nmin or nmax");
    }

    let ndata: usize = nmax - nmin;

    for i in (nmin..nmax-2).step_by(2) {
         s1 += func[i];
         s0 += func[i+1];
         s2 += func[i + 2];
     }

     //println!("sm1 =  {}, s0 = {}, sp1 = {}", s1, s0, s2);
     let mut res = step * (s1 + 4.0 * s0 + s2) / 3.0;
     //println!("step =  {}, res = {}", step, res);

     // If n is even, add the last slice separately
     if ndata % 2 == 0 {
        res += step * (5.0 * func[nmax-1] + 8.0 * func[nmax - 2] - func[nmax - 3]) / 12.0;
     }

    return res;

}


fn main() {
    const C2: f64 = 1.0/627.51; 
    const C9: f64 = 1.0e8/0.5291772; 
    const RGAS: f64 = 8.3144598/1000.0/2625.5; 
    const CLIGHT: f64 = 137.035999074; 

    let freq =  568.0; //cm-1
    let vf = 19.00*C2;
    let vb = 40.0*C2;
    let v0 = vf;

    let temp =  200.0; //K
    let beta = 1.0 / (RGAS*temp);
    let de = 0.001*C2; // from kcal/mol to Hartee
    let emax = 1000.0*C2; // from kcal/mol to Hartree;
    let omega = freq*CLIGHT/C9*2.0*PI;




    let kappa_wigner = wigner(beta, omega);

    let kappa_bell = bell(beta, omega);

    let kappa_st = skodje_truhlar(beta, omega, v0);

    let kappa_st_exact = skodje_truhlar_exact(beta, omega, v0);

    println!("Tunneling corrections at : {:>12.1} K", temp);
    println!("kappa Wigner    : {:>12.4}", kappa_wigner);
    println!("kappa Bell      : {:>12.4}", kappa_bell);
    println!("kappa Skodje    : {:>12.4}", kappa_st);
    println!("kappa Skodje ex : {:>12.4}", kappa_st_exact);

    let (kappa1, kappa2) = eckart(beta, omega, vf, vb, de, emax);

    println!("kappa Eckart 1  : {:>12.4}", kappa1);
    println!("kappa Eckart 2  : {:>12.4}", kappa2);
}
