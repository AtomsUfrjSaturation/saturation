
mod root;

use root::{find_interval, find_root};
use std::f64::consts::{E, PI};
use ndarray::{Array2, Array3, Array4, s};

const NAVO: f64 = 6.02214076e23; // Avogadro constant
const KB: f64 = 1.380649e-23; // Boltzmann (J/K)


pub struct PcSaftEos {
    m1: Array2<f64>,
    sigma: Array2<f64>,
    epsilon_k: Array2<f64>,
    m2: Option<Array2<f64>>,
    kbi: Array2<f64>,
    kab_k: Option<Array2<f64>>,
    eab_k: Option<Array2<f64>>,
    s: Option<Array2<f64>>,
    ap: Array2<f64>,
    bp: Array2<f64>,
    ncomp: usize,
    deltasimplified: bool,
}

impl PcSaftEos {

    pub fn new(m1: Array2<f64>, sigma: Array2<f64>, epsilon_k: Array2<f64>, m2: Option<Array2<f64>>, kbi: Array2<f64>, kab_k: Option<Array2<f64>>, eab_k: Option<Array2<f64>>, s: Option<Array2<f64>>, deltasimplified: bool) -> PcSaftEos {
        let mut ap = Array2::zeros((7, 3));
        let mut bp = Array2::zeros((7, 3));
        let ncomp = m1.len();

        ap[[0, 0]] = 0.91056314451539e0;
        ap[[0, 1]] = -0.30840169182720e0;
        ap[[0, 2]] = -0.09061483509767e0;
        ap[[1, 0]] = 0.63612814494991e0;
        ap[[1, 1]] = 0.18605311591713e0;
        ap[[1, 2]] = 0.45278428063920e0;
        ap[[2, 0]] = 2.68613478913903e0;
        ap[[2, 1]] = -2.50300472586548e0;
        ap[[2, 2]] = 0.59627007280101e0;
        ap[[3, 0]] = -26.5473624914884e0;
        ap[[3, 1]] = 21.4197936296668e0;
        ap[[3, 2]] = -1.72418291311787e0;
        ap[[4, 0]] = 97.7592087835073e0;
        ap[[4, 1]] = -65.2558853303492e0;
        ap[[4, 2]] = -4.13021125311661e0;
        ap[[5, 0]] = -159.591540865600e0;
        ap[[5, 1]] = 83.3186804808856e0;
        ap[[5, 2]] = 13.7766318697211e0;
        ap[[6, 0]] = 91.2977740839123e0;
        ap[[6, 1]] = -33.7469229297323e0;
        ap[[6, 2]] = -8.67284703679646e0;

        bp[[0, 0]] = 0.72409469413165e0;
        bp[[0, 1]] = -0.57554980753450e0;
        bp[[0, 2]] = 0.09768831158356e0;
        bp[[1, 0]] = 1.11913959304690e0 * 2e0;
        bp[[1, 1]] = 0.34975477607218e0 * 2e0;
        bp[[1, 2]] = -0.12787874908050e0 * 2e0;
        bp[[2, 0]] = -1.33419498282114e0 * 3e0;
        bp[[2, 1]] = 1.29752244631769e0 * 3e0;
        bp[[2, 2]] = -3.05195205099107e0 * 3e0;
        bp[[3, 0]] = -5.25089420371162e0 * 4e0;
        bp[[3, 1]] = -4.30386791194303e0 * 4e0;
        bp[[3, 2]] = 5.16051899359931e0 * 4e0;
        bp[[4, 0]] = 5.37112827253230e0 * 5e0;
        bp[[4, 1]] = 38.5344528930499e0 * 5e0;
        bp[[4, 2]] = -7.76088601041257e0 * 5e0;
        bp[[5, 0]] = 34.4252230677698e0 * 6e0;
        bp[[5, 1]] = -26.9710769414608e0 * 6e0;
        bp[[5, 2]] = 15.6044623461691e0 * 6e0;
        bp[[6, 0]] = -50.8003365888685e0 * 7e0;
        bp[[6, 1]] = -23.6010990650801e0 * 7e0;
        bp[[6, 2]] = -4.23812936930675e0 * 7e0;

        PcSaftEos {
            m1,
            sigma,
            epsilon_k,
            m2,
            kbi,
            kab_k,
            eab_k,
            s,
            ncomp,
            ap,
            bp,
            deltasimplified,
        }
    }

    pub fn pc_saft_d_t(&self, t: f64) -> Array2<f64> {
        let mut d_t = Array2::zeros((self.ncomp, 1));
        let sigma: &Array2<f64> = &self.sigma;
        let epsilon_k: &Array2<f64> = &self.epsilon_k;

        for i in 0..self.ncomp {
            d_t[[i, 0]] = sigma[[i, 0]] * (1.0 - 0.12 * (-3.0 * epsilon_k[[i, 0]] / t).exp());
        }
        d_t
    }

    pub fn pc_saft_csi(&self, dens: f64, t: f64, x: &Array2<f64>) -> Array2<f64> {
        let d_t = &self.pc_saft_d_t(t);
        let rho = &self.pc_saft_rho(dens);

        let mut csi: Array2<f64> = Array2::zeros((4, 1));

        for i in 0..4 {
            let mut sum = 0.0;
            for j in 0..self.ncomp {
                sum += x[[j, 0]] * self.m1[[j, 0]] * d_t[[j, 0]].powf(i as f64);
            }

            csi[[i, 0]] = sum * PI * rho / 6.0;
        }
        csi
    }

    pub fn pc_saft_ghs(&self, dens: f64, t: f64, x: &Array2<f64>) -> Array2<f64> {
        let csi = &self.pc_saft_csi(dens, t, x);
        let d_t = &self.pc_saft_d_t(t);

        let mut ghs: Array2<f64> = Array2::zeros((self.ncomp, self.ncomp));

        for i in 0..self.ncomp {
            for j in 0..self.ncomp {
                ghs[[i, j]] = 1.0 / (1.0 - csi[[3, 0]]) + (d_t[[i, 0]] * d_t[[j, 0]] / (d_t[[i, 0]] + d_t[[j, 0]])) * 3.0 * csi[[2, 0]] / ((1.0 - csi[[3, 0]]).powi(2)) + ((d_t[[i, 0]] * d_t[[j, 0]] / (d_t[[i, 0]] + d_t[[j, 0]])).powi(2)) * 2.0 * csi[[2, 0]].powi(2) / ((1.0 - csi[[3, 0]]).powi(3));
            }
        }
        ghs
    }

    pub fn pc_saft_zhs(&self, dens: f64, t: f64, x: &Array2<f64>) -> f64 {
        let csi = &self.pc_saft_csi(dens, t, x);

        let p1 = csi[[3, 0]] / (1.0 - csi[[3, 0]]);
        let p2 = 3.0 * csi[[1, 0]] * csi[[2, 0]] / (csi[[0, 0]] * (1.0 - csi[[3, 0]]).powi(2));
        let p3 = (3.0 * csi[[2, 0]].powi(3) - csi[[3, 0]] * csi[[2, 0]].powi(3)) / (csi[[0, 0]] * (1.0 - csi[[3, 0]]).powi(3));

        let zhs = p1 + p2 + p3;
        zhs
    }

    pub fn pc_saft_mmed(&self, x: &Array2<f64>) -> f64 {
        let mut mmed: f64 = 0.0;
        for i in 0..self.ncomp {
            mmed += x[[i, 0]] * self.m1[[i, 0]];
        }
        mmed
    }

    pub fn pc_saft_rho_dghsd_drho(&self, dens: f64, t: f64, x: &Array2<f64>) -> Array2<f64> {
        let csi = &self.pc_saft_csi(dens, t, x);
        let d_t = &self.pc_saft_d_t(t);

        let mut rho_dghsd_drho: Array2<f64> = Array2::zeros((self.ncomp, self.ncomp));

        for i in 0..self.ncomp {
            for j in 0..self.ncomp {
                rho_dghsd_drho[[i, j]] = csi[[3, 0]] / (1.0 - csi[[3, 0]]).powi(2) + (d_t[[i, 0]] * d_t[[j, 0]] / (d_t[[i, 0]] + d_t[[j, 0]])) * (3.0 * csi[[2, 0]] / (1.0 - csi[[3, 0]]).powi(2) + 6.0 * csi[[2, 0]] * csi[[3, 0]] / (1.0 - csi[[3, 0]]).powi(3)) + (d_t[[i, 0]] * d_t[[j, 0]] / (d_t[[i, 0]] + d_t[[j, 0]])).powi(2) * (4.0 * csi[[2, 0]].powi(2) / (1.0 - csi[[3, 0]]).powi(3) + 6.0 * csi[[2, 0]].powi(2) * csi[[3, 0]] / (1.0 - csi[[3, 0]]).powi(4));
            }
        }
        rho_dghsd_drho
    }

    pub fn pc_saft_zhc(&self, dens: f64, t: f64, x: &Array2<f64>) -> f64 {
        let mmed = &self.pc_saft_mmed(x);
        let zhs = &self.pc_saft_zhs(dens, t, x);
        let ghs = &self.pc_saft_ghs(dens, t, x);
        let rho_dghsd_drho = &self.pc_saft_rho_dghsd_drho(dens, t, x);

        let mut sum = 0.0;

        for i in 0..self.ncomp {
            sum += x[[i, 0]] * (self.m1[[i, 0]] - 1.0) / ghs[[i, i]] * rho_dghsd_drho[[i, i]];
        }
        let zhc = mmed * zhs - sum;
        zhc
    }

    pub fn pc_saft_a_e_b(&self, x: &Array2<f64>) -> (Array2<f64>, Array2<f64>) {
        let ap = &self.ap;
        let bp = &self.bp;
        let mmed = &self.pc_saft_mmed(x);

        let mut a: Array2<f64> = Array2::zeros((7, 1));
        let mut b: Array2<f64> = Array2::zeros((7, 1));

        for i in 0..7 {
            a[[i, 0]] = ap[[i, 0]] + (mmed - 1.0) * ap[[i, 1]] / mmed + (1.0 - 1.0 / mmed) * (1.0 - 2.0 / mmed) * ap[[i, 2]];
            b[[i, 0]] = bp[[i, 0]] + (mmed - 1.0) * bp[[i, 1]] / mmed + (1.0 - 1.0 / mmed) * (1.0 - 2.0 / mmed) * bp[[i, 2]];
        }
        (a, b)
    }

    pub fn pc_saft_c1(&self, dens: f64, t: f64, x: &Array2<f64>) -> f64 {
        let mmed = &self.pc_saft_mmed(x);
        let eta = &self.pc_saft_csi(dens, t, x)[[3, 0]];

        let c1 = (1.0 + mmed * (8.0 * eta - 2.0 * eta.powi(2)) / (1.0 - eta).powi(4) + (1.0 - mmed) * (20.0 * eta - 27.0 * eta.powi(2) + 12.0 * eta.powi(3) - 2.0 * eta.powi(4)) / ((1.0 - eta) * (2.0 - eta)).powi(2)).powi(-1);
        c1
    }

    pub fn pc_saft_c2(&self, dens: f64, t: f64, x: &Array2<f64>) -> f64 {
        let c1 = &self.pc_saft_c1(dens, t, x);
        let mmed = &self.pc_saft_mmed(x);
        let eta = &self.pc_saft_csi(dens, t, x)[[3, 0]];

        let c2 = -c1.powi(2) * (mmed * (-4.0 * eta.powi(2) + 20.0 * eta + 8.0) / (1.0 - eta).powi(5) + (1.0 - mmed) * (2.0 * eta.powi(3) + 12.0 * eta.powi(2) - 48.0 * eta + 40.0) / ((1.0 - eta) * (2.0 - eta)).powi(3));
        c2
    }

    pub fn pc_saft_i1_e_i2(&self, dens: f64, t: f64, x: &Array2<f64>) -> (f64, f64) {
        let (a, b) = &self.pc_saft_a_e_b(x);
        let eta = &self.pc_saft_csi(dens, t, x)[[3, 0]];

        let mut i1: f64 = 0.0;
        let mut i2: f64 = 0.0;

        for i in 0..7 {
            i1 += a[[i, 0]] * eta.powf(i as f64);
            i2 += b[[i, 0]] * eta.powf(i as f64);
        }

        (i1, i2)
    }

    pub fn pc_saft_detai1_deta_e_detai2_deta(&self, dens: f64, t: f64, x: &Array2<f64>) -> (f64, f64) {
        let (a, b) = &self.pc_saft_a_e_b(x);
        let eta = &self.pc_saft_csi(dens, t, x)[[3, 0]];

        let mut deta_i1_deta: f64 = 0.0;
        let mut deta_i2_deta: f64 = 0.0;

        for i in 0..7 {
            deta_i1_deta += a[[i, 0]] * (i as f64 + 1.0) * eta.powf(i as f64);
            deta_i2_deta += b[[i, 0]] * (i as f64 + 1.0) * eta.powf(i as f64);
        }

        (deta_i1_deta, deta_i2_deta)
    }

    pub fn pc_saft_mat_sigma(&self) -> Array2<f64> {
        let sigma = &self.sigma;
        let mut mat_sigma: Array2<f64> = Array2::zeros((self.ncomp, self.ncomp));

        for i in 0..self.ncomp {
            for j in 0..self.ncomp {
                mat_sigma[[i, j]] = (sigma[[i, 0]] + sigma[[j, 0]]) / 2.0;
            }
        }
        mat_sigma
    }

    pub fn pc_saft_mat_epsilon_k(&self, x: &Array2<f64>) -> Array2<f64> {
        let epsilon_k = &self.epsilon_k;
        let kbi = &self.kbi;
        let mut mat_epsilon_k = Array2::zeros((self.ncomp, self.ncomp));

        for i in 0..self.ncomp {
            for j in 0..self.ncomp {
                mat_epsilon_k[[i, j]] = (epsilon_k[[j, 0]] * epsilon_k[[j, 0]]).powf(0.5) * (1.0 - kbi[[i, j]]);
            }
        }
        mat_epsilon_k
    }

    pub fn pc_saft_m_2esig_3_e_m_2e_2sig_3(&self, t: f64, x: &Array2<f64>) -> (f64, f64) {
        let m1 = &self.m1;
        let mat_epsilon_k = self.pc_saft_mat_epsilon_k(x);
        let mat_sigma = self.pc_saft_mat_sigma();
        let mut m_2esig_3: f64 = 0.0;
        let mut m_2e_2sig_3: f64 = 0.0;

        for i in 0..self.ncomp {
            for j in 0..self.ncomp {
                m_2esig_3 += x[[i, 0]] * x[[j, 0]] * m1[[i, 0]] * m1[[j, 0]] * mat_epsilon_k[[i, j]] / t * (mat_sigma[[i, j]]).powi(3);
                m_2e_2sig_3 += x[[i, 0]] * x[[j, 0]] * m1[[i, 0]] * m1[[j, 0]] * (mat_epsilon_k[[i, j]] / t).powi(2) * (mat_sigma[[i, j]]).powi(3);
            }
        }
        (m_2esig_3, m_2e_2sig_3)
    }

    pub fn pc_saft_zdisp(&self, dens: f64, t: f64, x: &Array2<f64>) -> f64 {
        let mmed = &self.pc_saft_mmed(x);
        let (detai1_deta, detai2_deta) = &self.pc_saft_detai1_deta_e_detai2_deta(dens, t, x);
        let c1 = &self.pc_saft_c1(dens, t, x);
        let c2 = &self.pc_saft_c2(dens, t, x);
        let eta = &self.pc_saft_csi(dens, t, x)[[3, 0]];
        let (i1, i2) = &self.pc_saft_i1_e_i2(dens, t, x);
        let (m_2esig_3, m_2e_2sig_3) = &self.pc_saft_m_2esig_3_e_m_2e_2sig_3(t, x);
        let rho = &self.pc_saft_rho(dens);

        let zdisp = -2.0 * PI * rho * detai1_deta * m_2esig_3 - PI * rho * mmed * (c1 * detai2_deta + c2 * eta * i2) * m_2e_2sig_3;
        zdisp
    }

    pub fn pc_saft_z(&self, dens: f64, t: f64, x: &Array2<f64>) -> f64 {
        let zhc = &self.pc_saft_zhc(dens, t, x);
        let zdisp = &self.pc_saft_zdisp(dens, t, x);
        let s = &self.s;

        match s {
            Some(k) => {
                let da_ass_deta = self.pc_saft_da_ass_deta(dens, t, x);
                let zassoc = dens * da_ass_deta;
                let z = 1. + zhc + zdisp + zassoc;
                z
            },
            None => {
                let z = 1.0 + zhc + zdisp;
                z
            }
        }
    }

    pub fn pc_saft_pressure(&self, dens: f64, t: f64, x: &Array2<f64>) -> f64 {
        let z = self.pc_saft_z(dens, t, x);
        let p = z * KB * t * dens * NAVO;
        p
    }

    pub fn pc_saft_rho(&self, dens: f64) -> f64 {
        let rho = dens * NAVO / 1.0e30;
        rho
    }

    pub fn pc_saft_dens(&self, t: f64, p: f64, x: &Array2<f64>, phase: Option<&str>) -> (f64, f64) {
        let m1 = &self.m1;
        let d_t = self.pc_saft_d_t(t);
        let mut soma = 0.;

        // for liquid
        let etaguessl = 0.5;
        // for gas
        let etaguessv = 1e-10;
        
        for i in 0..self.ncomp {
            soma += x[[i, 0]] * m1[[i, 0]] * d_t[[i, 0]].powi(3);
        }
        
        // Intial guesses of the optimization
        let densl0 = 6.0 / PI * etaguessl / soma * 1e30 / NAVO;
        let densv0 = 6.0 / PI * etaguessv / soma * 1e30 / NAVO;
        
        // Now, finding the root for liquid and vapor:
        let residuo_liq = |dens: f64| -> f64 {
            let f = (self.pc_saft_pressure(dens, t, x)/p)-1.;
            f
        };


        let residuo_vap = |dens: f64| -> f64 {
            let f = f64::log(self.pc_saft_pressure(dens, t, x)/p, E);
            f
        };
            
        
        let interval_liq = find_interval(densl0, residuo_liq);
        let interval_vap = find_interval(densv0, residuo_vap);
        
        match phase {
            Some(k) => {
                if k == "liq" {
                    (find_root(interval_liq, residuo_liq), 0.)
                } else if k == "vap" {
                    (0., find_root(interval_vap, residuo_vap))
                } else {
                    println!("{}", "Invalid Phase!!");
                    (0., 0.)
                }
            },
            None => {
                (find_root(interval_liq, residuo_liq), find_root(interval_vap, residuo_vap))
            }
        }
    }

    pub fn pc_saft_a_hs(&self, dens: f64, t: f64, x: &Array2<f64>) -> f64 {
        let csi = self.pc_saft_csi(dens, t, x);
        let a_hs = (3. * csi[[1, 0]] * csi[[2, 0]] / (1. - csi[[3, 0]])
            + csi[[2, 0]].powi(3) / (csi[[3, 0]] * (1. - csi[[3, 0]]).powi(2))
            + (csi[[2, 0]].powi(3) / csi[[3, 0]].powi(2) - csi[[0, 0]]) * f64::log(1. - csi[[3, 0]], E)) / csi[[0, 0]];
        a_hs
    }

    pub fn pc_saft_a_hc(&self, dens: f64, t: f64, x: &Array2<f64>) -> f64 {
        let mmed = self.pc_saft_mmed(x);
        let m1 = &self.m1;
        let ghs = self.pc_saft_ghs(dens, t, x);
        let a_hs = self.pc_saft_a_hs(dens, t, x);
        let mut soma = 0.;

        for i in 0..self.ncomp {
            soma += -x[[i, 0]] * (m1[[i, 0]] - 1.) * f64::log(ghs[[i, 0]], E);
        }

        let a_hc = mmed * a_hs + soma;
        a_hc
    }

    pub fn pc_saft_a_disp(&self, dens: f64, t: f64, x: &Array2<f64>) -> f64 {
        let (i1, i2) = self.pc_saft_i1_e_i2(dens, t, x);
        let (m_2esig_3, m_2e_2sig_3) = self.pc_saft_m_2esig_3_e_m_2e_2sig_3(t, x);
        let rho = self.pc_saft_rho(dens);
        let mmed = self.pc_saft_mmed(x);
        let c1 = self.pc_saft_c1(dens, t, x);

        let a_disp = -2. * PI * rho * i1 * m_2esig_3 -
            PI * rho * mmed * c1 * i2 * m_2e_2sig_3;
        a_disp
    }

    pub fn pc_saft_a_res(&self, dens: f64, t: f64, x: &Array2<f64>) -> f64 {
        match &self.s {
            Some(_) => {
                let a_ass = self.pc_saft_a_ass(dens, t, x).unwrap_or(0.);
                let a_res = self.pc_saft_a_hc(dens, t, x) + self.pc_saft_a_disp(dens, t, x) + a_ass;
                a_res
            }
            None => {
                let a_res = self.pc_saft_a_hc(dens, t, x) + self.pc_saft_a_disp(dens, t, x);
                a_res
            }
        }
    }

    pub fn pc_saft_mat_dcsi_dxk(&self, dens: f64, t: f64) -> Array2<f64> {
        let m1 = &self.m1;
        let d_t = self.pc_saft_d_t(t);
        let rho = self.pc_saft_rho(dens);
        let mut mat_dcsi_dxk = Array2::<f64>::zeros((4, self.ncomp));


        for k in 0..self.ncomp {
            for i in 0..4 {
                mat_dcsi_dxk[[i, k]] = PI / 6. * rho * m1[[k, 0]] * d_t[[k, 0]].powi(i as i32);
            }
        }
        mat_dcsi_dxk
    }

    pub fn pc_saft_dghs_dx(&self, dens: f64, t: f64, x: &Array2<f64>) -> Array2<Array2<f64>> {
        let mat_dcsi_dxk = self.pc_saft_mat_dcsi_dxk(dens, t);
        let csi = self.pc_saft_csi(dens, t, x);
        let d_t = self.pc_saft_d_t(t);

        let mut vector = vec!();
        let y: Array2<f64> = Array2::zeros((self.ncomp, self.ncomp));

        for i in 0..self.ncomp {
            vector.push(y.clone());
        }

        let mut dghs_dx: Array2<Array2<f64>> = Array2::from_shape_vec((self.ncomp, 1), vector).unwrap();

        for k in 0..self.ncomp {
            let p1 = mat_dcsi_dxk[[3, k]] / (1. - csi[[3, 0]]).powi(2);
            let p2 = 3. * mat_dcsi_dxk[[2, k]] / (1. - csi[[3, 0]]).powi(2) +
                6. * (csi[[2, 0]] * mat_dcsi_dxk[[3, k]]) / (1. - csi[[3, 0]]).powi(3);
            let p3 = 4. * csi[[2, 0]] * mat_dcsi_dxk[[2, k]] / (1. - csi[[3, 0]]).powi(3) +
                6. * csi[[2, 0]].powi(2) * mat_dcsi_dxk[[3, k]] / (1. - csi[[3, 0]]).powi(4);
            for i in 0..self.ncomp {
                for j in 0..self.ncomp {
                    dghs_dx[[k, 0]][[i, j]] = p1 + (d_t[[i, 0]] * d_t[[j, 0]] / (d_t[[i, 0]] + d_t[[j, 0]])) * p2 +
                        (d_t[[i, 0]] * d_t[[j, 0]] / (d_t[[i, 0]] + d_t[[j, 0]])).powi(2) * p3;
                }
            }
        }
        dghs_dx
    }

    pub fn pc_saft_da_hs_dx(&self, dens: f64, t: f64, x: &Array2<f64>) -> Array2<f64> {
        let mat_dcsi_dxk = self.pc_saft_mat_dcsi_dxk(dens, t);
        let a_hs = self.pc_saft_a_hs(dens, t, x);
        let csi = self.pc_saft_csi(dens, t, x);
        let mut da_hs_dx: Array2<f64> = Array2::zeros((self.ncomp, 1));

        for i in 0..self.ncomp {
            let p1 = -mat_dcsi_dxk[[0, i]] / csi[[0, 0]] * a_hs;
            let p2 = 3. * (mat_dcsi_dxk[[1, i]] * csi[[2, 0]]
                + csi[[1, 0]] * mat_dcsi_dxk[[2, i]]) / (1. - csi[[3, 0]]);
            let p3 = 3. * csi[[1, 0]] * csi[[2, 0]] * mat_dcsi_dxk[[3, i]
                ] / (1. - csi[[3, 0]]).powi(2);
            let p4 = 3. * csi[[2, 0]].powi(2) * mat_dcsi_dxk[[2, i]]
                / (csi[[3, 0]] * (1. - csi[[3, 0]]).powi(2));
            let p5 = csi[[2, 0]].powi(3) * mat_dcsi_dxk[[3, i]]
                * (3. * csi[[3, 0]] - 1.) / (csi[[3, 0]].powi(2)
                * (1. - csi[[3, 0]]).powi(3));
            let p6 = ((3. * csi[[2, 0]].powi(2) * mat_dcsi_dxk[[2, i]]
                * csi[[3, 0]] - 2. * csi[[2, 0]].powi(3) * mat_dcsi_dxk[[3, i]])
                / csi[[3, 0]].powi(3) - mat_dcsi_dxk[[0, i]])
                * (1. - csi[[3, 0]]).ln();
            let p7 = (csi[[0, 0]] - csi[[2, 0]].powi(3)
                / csi[[3, 0]].powi(2)) * mat_dcsi_dxk[[3, i]]
                / (1. - csi[[3, 0]]);
            da_hs_dx[[i, 0]] = p1 + (p2 + p3 + p4 + p5 + p6 + p7) / csi[[0, 0]]
        }
        da_hs_dx
    }

    pub fn pc_saft_da_hc_dx(&self, dens: f64, t: f64, x: &Array2<f64>) -> Array2<f64> {
        let m1 = &self.m1;
        let mmed = self.pc_saft_mmed(x);
        let a_hs = self.pc_saft_a_hs(dens, t, x);
        let ghs = self.pc_saft_ghs(dens, t, x);
        let dghs_dx = self.pc_saft_dghs_dx(dens, t, x);
        let da_hs_dx = self.pc_saft_da_hs_dx(dens, t, x);
        let mut da_hc_dx: Array2<f64> = Array2::zeros((self.ncomp, 1));
        let mut somai: Array2<f64> = Array2::zeros((self.ncomp, 1));

        for k in 0..self.ncomp {
            for i in 0..self.ncomp {
                somai[[k, 0]] = x[[i, 0]] * (m1[[i, 0]] - 1.) / ghs[[i, i]] * dghs_dx[[k, 0]][[i, i]];
            }
            da_hc_dx[[k, 0]] = m1[[k, 0]] * a_hs + mmed * da_hs_dx[[k, 0]] - somai[[k, 0]] 
            - (m1[[k, 0]] - 1.) * f64::log(ghs[[k, k]], E);
        }
        da_hc_dx
    }

    pub fn pc_saft_dai_dx_e_dbi_dx(&self, x: &Array2<f64>) -> (Array2<Array2<f64>>, Array2<Array2<f64>>) {
        let m1 = &self.m1;
        let mmed = self.pc_saft_mmed(x);
        let ap = &self.ap;
        let bp = &self.bp;
        let mut dai_dx: Array2<f64>;
        let mut dbi_dx: Array2<f64>;
        let mut vector = Vec::new();

        let y: Array2<f64> = Array2::zeros((7, 1));

        for i in 0..self.ncomp {
            vector.push(y.clone());
        }

        let mut dai_dx: Array2<Array2<f64>> = Array2::from_shape_vec((self.ncomp, 1), vector.clone()).unwrap();
        let mut dbi_dx: Array2<Array2<f64>> = Array2::from_shape_vec((self.ncomp, 1), vector.clone()).unwrap();

        for k in 0..self.ncomp {
            for i in 0..7 {
                dai_dx[[k, 0]][[i, 0]] = (m1[[k, 0]] / mmed.powi(2)) *
                    ap[[i, 1]] + m1[[k, 0]] / (mmed.powi(2)) * (3. - 4. / mmed) * ap[[i, 2]];
                dbi_dx[[k, 0]][[i, 0]] = (m1[[k, 0]] / mmed.powi(2)) *
                    bp[[i, 1]] + m1[[k, 0]] / (mmed.powi(2)) * (3. - 4. / mmed) * bp[[i, 2]];
            }
        }
        (dai_dx, dbi_dx)
    }

    pub fn pc_saft_di1_dx_e_di2_dx(&self, dens: f64, t: f64, x: &Array2<f64>) -> (Array2<f64>, Array2<f64>) {
        let (a, b) = self.pc_saft_a_e_b(x);
        let eta = self.pc_saft_csi(dens, t, x)[[3, 0]];
        let (dai_dx, dbi_dx) = self.pc_saft_dai_dx_e_dbi_dx(x);
        let mat_dcsi_dxk = self.pc_saft_mat_dcsi_dxk(dens, t);
        let mut di1_dx = Array2::zeros((self.ncomp, 1));
        let mut di2_dx = Array2::zeros((self.ncomp, 1));

        for k in 0..self.ncomp {
            for i in 0..7 {
                di1_dx[[k, 0]] += a[[i, 0]] * (i as f64) * mat_dcsi_dxk[[3, k]] * eta.powi(i as i32 - 1) + dai_dx[[k, 0]][[i, 0]] * eta.powi(i as i32);
                di2_dx[[k, 0]] += b[[i, 0]] * (i as f64) * mat_dcsi_dxk[[3, k]] * eta.powi(i as i32 - 1) + dbi_dx[[k, 0]][[i, 0]] * eta.powi(i as i32);
            }
        }
        (di1_dx, di2_dx)
    }

    pub fn pc_saft_dc1_dx(&self, dens: f64, t: f64, x: &Array2<f64>) -> Array2<f64> {
        let m1 = &self.m1;
        let eta = self.pc_saft_csi(dens, t, x)[[3, 0]];
        let mat_dcsi_dxk = self.pc_saft_mat_dcsi_dxk(dens, t);
        let c1 = self.pc_saft_c1(dens, t, x);
        let c2 = self.pc_saft_c2(dens, t, x);
        let mut dc1_dx = Array2::zeros((self.ncomp, 1));

        for k in 0..self.ncomp {
            let p1 = c2 * mat_dcsi_dxk[[3, k]];
            let p2 = m1[[k, 0]] * (8. * eta - 2. * eta.powi(2)) / (1. - eta).powi(4);
            let p3 = -m1[[k, 0]] * (20. * eta - 27. * eta.powi(2) + 12. * eta.powi(3) - 2.
                * eta.powi(4)) / ((1. - eta) * (2. - eta)).powi(2);
            dc1_dx[[k, 0]] = p1 - c1.powi(2) * (p2 + p3);
        }
        dc1_dx
    }

    pub fn pc_saft_dm_2esig_3_dx_e_dm_2e_2sig_3_dx(&self, dens: f64, t: f64, x: &Array2<f64>) -> (Array2<f64>, Array2<f64>){
        let m1 = &self.m1;
        let mat_epsilon_k = self.pc_saft_mat_epsilon_k(x);
        let mat_sigma = self.pc_saft_mat_sigma();
        let mut dm_2esig_3_dx:Array2<f64> = Array2::zeros((self.ncomp, 1));
        let mut dm_2e_2sig_3_dx:Array2<f64> = Array2::zeros((self.ncomp, 1));

        for k in 0..self.ncomp{
            for i in 0..self.ncomp{
                dm_2esig_3_dx[[k,0]] += 2.0*m1[[k,0]]*x[[i,0]]*m1[[i,0]] * (mat_epsilon_k[[k, i]]/t) * (mat_sigma[[k, i]]).powi(3);
                dm_2e_2sig_3_dx[[k,0]] += 2.0*m1[[k,0]]*x[[i,0]]*m1[[i,0]] * (mat_epsilon_k[[k, i]]/t).powi(2) * (mat_sigma[[k, i]]).powi(3);
            }
        }

        return (dm_2esig_3_dx, dm_2e_2sig_3_dx);
    }

    pub fn pc_saft_da_disp_dx(&self, dens: f64, t: f64, x: &Array2<f64>) -> Array2<f64> {
        let m1 = &self.m1;
        let (dI1_dx, dI2_dx) = self.pc_saft_di1_dx_e_di2_dx(dens, t, x);
        let (dm_2esig_3_dx, dm_2e_2sig_3_dx) = self.pc_saft_dm_2esig_3_dx_e_dm_2e_2sig_3_dx(dens, t, x);
        let dC1_dx = self.pc_saft_dc1_dx(dens, t, x);
        let C1 = self.pc_saft_c1(dens, t, x);
        let (I1, I2) = self.pc_saft_i1_e_i2(dens, t, x);
        let (m_2esig_3, m_2e_2sig3) = self.pc_saft_m_2esig_3_e_m_2e_2sig_3(t, x);
        let mmed = self.pc_saft_mmed(x);
        let rho= self.pc_saft_rho(dens);
        let mut da_disp_dx: Array2<f64> = Array2::zeros((self.ncomp, 1));

        for k in 0..self.ncomp{
            let p1 = -2.*PI*rho*(dI1_dx[[k,0]]*m_2esig_3 + I1*dm_2esig_3_dx[[k,0]]);
            let p2 = (m1[[k,0]]*C1*I2 + mmed*dC1_dx[[k,0]]*I2+mmed*C1*dI2_dx[[k,0]])*m_2e_2sig3;
            let p3 = mmed*C1*I2*dm_2e_2sig_3_dx[[k,0]];
            da_disp_dx[[k,0]] = p1 - PI*rho*(p2+p3);
        }
        da_disp_dx
    }

    pub fn pc_saft_da_ass_dx_num(&self, dens: f64, t: f64, x: &Array2<f64>) -> Result<Array2<f64>, &'static str> {
        let step = 1e-4;
        let mut da_ass_dx_vec: Vec<f64> = vec![];

        match &self.s {
            Some(s_ref) => {
                let nsite = s_ref.len();
                let mut counter = 0;

                for k in 0..self.ncomp {
                    let mut xmais = x.clone();
                    let mut xmenos = x.clone();

                    xmais[[k, 0]]  = xmais[[k, 0]] + step;
                    xmenos[[k, 0]] = xmenos[[k, 0]] - step;

                    let x_a_mais = self.pc_saft_x_tan(dens, t, &xmais)?;
                    let x_a_menos = self.pc_saft_x_tan(dens, t, &xmenos)?;

                    let mut a_ass_mais  = 0.;
                    let mut a_ass_menos = 0.;

                    for i in 0..self.ncomp {
                        let mut somamais  = 0.;
                        let mut somamenos = 0.;

                        for j in 0..nsite {
                            somamais += f64::log(x_a_mais[[j, i]], E) - x_a_mais[[j, i]]/2.;
                            somamenos += f64::log(x_a_menos[[j, i]], E) - x_a_menos[[j, i]]/2.;
                        }

                        a_ass_mais += xmais[[i, 0]] * (somamais + s_ref.column(i).sum() * 0.5);
                        a_ass_menos += xmenos[[i,0]] * (somamenos + s_ref.column(i).sum() * 0.5);
                    }
                    counter += 1;
                    da_ass_dx_vec.push((a_ass_mais - a_ass_menos)/2./step);
                }

                let da_ass_dx: Array2<f64> = match Array2::from_shape_vec((counter, 1), da_ass_dx_vec) {
                    Ok(array) => array,
                    Err(e) => panic!("Falha na conversão: {}", e),
                };
                Ok(da_ass_dx)
            },
            None => Err("A molécula não possui associação e essa função não pode ser chamada para ela"),
        }
    
    }
    
    pub fn pc_saft_mu_assoc_kt(&self, dens: f64, t: f64, x: &Array2<f64>) -> Result<Array2<f64>, &'static str> {
        match &self.s {
            Some(s_ref) => {

                let x_a = self.pc_saft_x_tan(dens, t, x)?;
                let rho = self.pc_saft_rho(dens);
                let dx_a_drho_k = self.pc_saft_dx_a_dx_k_num(dens, t, x).unwrap()/rho;
                
                let mut sum1: Array2<f64> = Array2::zeros((x_a.ncols(),1));
                let mut sum3: Array2<f64> = Array2::<f64>::zeros((dx_a_drho_k.dim().2, 1));
    
                for i in 0..x_a.ncols() {
                    for j in 0..x_a.nrows() {
                        sum1[[i,0]] += x_a[[j, i]].ln() - x_a[[j, i]] * 0.5 + 0.5 * s_ref[[j, i]];
                    }
                }
                
                for m in 0..dx_a_drho_k.dim().2 {
                    for i in 0..dx_a_drho_k.dim().1 {
                        for j in 0..dx_a_drho_k.dim().0 {
                            sum3[[m, 0]] += rho * x[[i, 0]] * dx_a_drho_k[[j, i, m]] * (1.0 / x_a[[j, i]] - 0.5);
                        }
                    }
                }

                let mut mu_assoc_kt: Array2<f64> = Array2::zeros((sum1.len(), 1));
                for i in 0..sum1.len(){
                    mu_assoc_kt[[i, 0]] = sum1[[i,0]] + sum3[[i,0]];
                }

                Ok(mu_assoc_kt)
            },
            None => Err("A molécula não possui associação e essa função não pode ser chamada para ela"),
        }
    }

    pub fn pc_saft_da_res_dx(&self, dens: f64, t: f64, x: &Array2<f64>) -> Result<Array2<f64>, &'static str> {
        let da_hc_dx = self.pc_saft_da_hc_dx(dens, t, x);
        let da_disp_dx = self.pc_saft_da_disp_dx(dens, t, x);
    
        match &self.s {
            Some(_) => {
                let da_ass_dx = self.pc_saft_da_ass_dx_num(dens, t, x)?;
                let da_res_dx: Array2<f64> = da_hc_dx + da_disp_dx + da_ass_dx;
                return Ok(da_res_dx);
            },
            None => return Ok(da_hc_dx + da_disp_dx),
        }
    }

    pub fn pc_saft_mu_res_kt(&self, dens: f64 , t: f64, x: &Array2<f64>) -> Array2<f64> {
        let mut mu_res_kt: Array2<f64> = Array2::zeros((self.ncomp, 1));
        let mut mu_hc: Array2<f64> = Array2::zeros((self.ncomp, 1));
        let mut mu_disp: Array2<f64> = Array2::zeros((self.ncomp, 1));
    
        let mu_assoc_kt = match self.pc_saft_mu_assoc_kt(dens, t, x) {
            Ok(array) => array,
            Err(_) => panic!(),
        };
        let zhc = self.pc_saft_zhc(dens, t, x);
        let zdisp = self.pc_saft_zdisp(dens, t, x);
        let da_hc_dx = self.pc_saft_da_hc_dx(dens, t, x);
        let da_disp_dx = self.pc_saft_da_disp_dx(dens, t, x);
        let ares_hc = self.pc_saft_a_hc(dens, t, x);
        let ares_disp = self.pc_saft_a_disp(dens, t, x);
    
        let somahc: f64 = x.iter()
                        .zip(da_hc_dx.iter())
                        .map(|(x_i, da_hc_dx_i)| x_i * da_hc_dx_i)
                        .sum();

        let somadisp: f64 = x.iter()
                        .zip(da_disp_dx.iter())
                        .map(|(x_i, da_disp_dx_i)| x_i * da_disp_dx_i)
                        .sum();

        for i in 0..self.ncomp {
            mu_hc[[i, 0]] = ares_hc + zhc + da_hc_dx[[i, 0]] - somahc;
            mu_disp[[i, 0]] = ares_disp + zdisp + da_disp_dx[[i, 0]] - somadisp;
            mu_res_kt[[i, 0]] = mu_hc[[i, 0]] + mu_disp[[i, 0]] + mu_assoc_kt[[i, 0]];
        }

        mu_res_kt
    }

    pub fn pc_saft_phi(&self, dens: f64, t: f64, x: &Array2<f64>) -> Array2<f64> {
        let mu = self.pc_saft_mu_res_kt(dens, t, x);
        let z = self.pc_saft_z(dens, t, x);
        let mut lnphi: Array2<f64> = Array2::zeros((self.ncomp, 1));

        for i in 0..self.ncomp {
            lnphi[[i, 0]] = mu[[i, 0]] - f64::log(z, E);
        }

        let mut phi:Array2<f64> = Array2::zeros((lnphi.len(), 1));
        for i in 0..phi.len(){
            phi[[i, 0]] = lnphi[[i, 0]].exp();
        }
        phi
    }

    pub fn pc_saft_eaibj_k(&self) -> Array2<f64> {
        let mut eaibj_k: Array2<f64> = Array2::zeros((self.ncomp, self.ncomp));
        let eab_k = &self.eab_k;

        match eab_k {
            Some(k) => {
                for i in 0..self.ncomp {
                    for j in 0..self.ncomp {
                        eaibj_k[[i, j]] = (k[[i, 0]] + k[[j, 0]])/2.;
                    }
                }
                eaibj_k
            },
            None => {
                eaibj_k
            }
        }
    }

    pub fn pc_saft_kaibj_k(&self) -> Array2<f64> {
        let mut kaibj_k: Array2<f64> = Array2::zeros((self.ncomp, self.ncomp));
        let kab_k = &self.kab_k;
        let mat_sigma = self.pc_saft_mat_sigma();

        match kab_k {
            Some(k) => {
                for i in 0..self.ncomp {
                    for j in 0..self.ncomp {
                        kaibj_k[[i, j]] = (k[[i, 0]]*k[[j, 0]]).sqrt() * ((mat_sigma[[i, i]]*mat_sigma[[j, j]]).sqrt() / ((mat_sigma[[i, i]] + mat_sigma[[j, j]])/2.)).powi(3);
                    }
                }
                kaibj_k
            },
            None => {
                kaibj_k
            }
        }
    }

    pub fn pc_saft_delt(&self, dens: f64, t: f64, x: &Array2<f64>) -> Result<Array4<f64>, &'static str> {
        match &self.s {
            Some(s) => {
                let nsite = s.len();
                let mat_sigma = self.pc_saft_mat_sigma();
                let e_ai_bj_k = self.pc_saft_eaibj_k();
                let k_ai_bj_k = self.pc_saft_kaibj_k();
                let ghs = self.pc_saft_ghs(dens, t, x);
                let d_t = self.pc_saft_d_t(t); 
                let mut delt: Array4<f64> = Array4::zeros((nsite, self.ncomp, nsite, self.ncomp));
                
                if self.deltasimplified == false {
                    for i in 0..self.ncomp {
                        for j in 0..nsite {
                            for k in 0..self.ncomp {
                                for l in 0..nsite {
                                    if j != l {
                                        delt[[j, i, l, k]] = (d_t[[i, 0]] * d_t[[k, 0]] / (d_t[[i, 0]] + d_t[[k, 0]])).powi(3) * 
                                        ghs[[i, k]] * ((e_ai_bj_k[[i, k]] / t).exp() - 1.0) * k_ai_bj_k[[i, k]];
                                    }
                                }
                            }
                        }
                    }
                } else {
                    for i in 0..self.ncomp {
                        for j in 0..nsite {
                            for k in 0..self.ncomp {
                                for l in 0..nsite {
                                    if j != l {
                                        delt[[j, i, l, k]] = mat_sigma[[k, i]].powi(3) * 
                                        ghs[[k, i]] * ((e_ai_bj_k[[k, i]] / t).exp() - 1.0) * k_ai_bj_k[[k, i]];
                                    }
                                }
                            }
                        }
                    }
                }
                Ok(delt)
            },
            None => Err("A molécula não possui associação e essa função não pode ser chamada para ela"),
        }
    }
    
    pub fn pc_saft_x_tan(&self, dens: f64, t: f64, x: &Array2<f64>) -> Result<Array2<f64>, &'static str> {
        match &self.s {
            Some(s_ref) => {
                let delta = self.pc_saft_delt(dens, t, x)?;
                let s = s_ref.clone();
                let nsite = s.len();
                let ncomp = self.ncomp;
                let rho = self.pc_saft_rho(dens);
                let mut x_a = Array2::ones((nsite, ncomp)) * 0.5;
                let mut res = 1.0;
                let mut it = 0;
                let it_max = 3000;
    
                while res > 1e-9 && it < it_max {
                    it += 1;
                    let x_a_old = x_a.clone();
                    let mut sum1: Array2<f64> = Array2::zeros((nsite, ncomp));
                    
                    for i in 0..nsite {
                        for j in 0..ncomp {
                            for k in 0..nsite {
                                for l in 0..ncomp {
                                    sum1[[i, j]] += s[[k, l]] * x_a_old[[k, l]] * delta[[k, l, i, j]] * x[[l, 0]];
                                }
                            }
                        }
                    }
    
                    x_a = 1.0 / (1.0 + rho * sum1);

                    let dif = &x_a - &x_a_old;
                    let mut res_temp = f64::NEG_INFINITY;
                    for &item in dif.iter() {
                        if item.abs() > res_temp {
                            res_temp = item.abs();
                        }
                    }
                    res = res_temp;
                }

                if it == it_max {
                    x_a *= f64::NAN;
                }
    
                Ok(x_a)
            },
            None => Err("A variável 's' não foi definida para esta molécula."),
        }
    }
    
    pub fn pc_saft_dx_a_dx_k_num(&self, dens: f64, t: f64, x: &Array2<f64>) -> Result<Array3<f64>, &'static str> {
        match &self.s {
            Some(s) => {
                let nsite = s.len();
                let step = 1e-8;
                let mut dx_a_dx_k = Array3::ones((nsite, self.ncomp, self.ncomp));
    
                for k in 0..self.ncomp {
                    let mut xmais = x.clone();
                    let mut xmenos = x.clone();
    
                    xmais[[k, 0]] += step;
                    xmenos[[k, 0]] -= step;
    
                    let x_a_mais = self.pc_saft_x_tan(dens, t, &xmais)?;
                    let x_a_menos = self.pc_saft_x_tan(dens, t, &xmenos)?;
    
                    for i in 0..nsite {
                        for j in 0..self.ncomp {
                            dx_a_dx_k[[i, j, k]] = (x_a_mais[[i, j]] - x_a_menos[[i, j]]) / (2.0 * step);
                        }
                    }
                }
    
                Ok(dx_a_dx_k)
            },
            None => Err("A variável 's' não foi definida para esta molécula."),
        }
    }

    pub fn pc_saft_a_ass(&self, dens: f64, t: f64, x: &Array2<f64>) -> Result<f64, &'static str> {
        match &self.s {
            Some(s_ref) => {
                let s = s_ref.clone();
                let nsite = s_ref.clone().len();
                let x_a = self.pc_saft_x_tan(dens, t, x)?;
                let mut a_ass = 0.;
    
                for i in 0..self.ncomp {
                    let mut s1 = 0.0;
                    for j in 0..nsite {
                        s1 += (x_a[[j, i]].ln() - x_a[[j, i]] / 2.0 + 0.5) * s_ref[[j,i]];
                    }
                    let s_slice = s.slice(s![.., i]);
                    let sum_s = s_slice.sum();
                    a_ass += x[[i, 0]] * (s1);
                }
    
                Ok(a_ass)
            },
            None => Err("A variável 's' não foi definida para esta molécula."),
        }
    }

    pub fn pc_saft_da_ass_deta(&self, dens: f64, t: f64, x: &Array2<f64>) -> f64 {
        let h = 1e-4;
        let a_ass_mais = self.pc_saft_a_ass(dens + h, t, x).unwrap_or(0.);
        let a_ass_menos = self.pc_saft_a_ass(dens - h, t, x).unwrap_or(0.);
        let res = (a_ass_mais - a_ass_menos) / (2. * h);
        res
    }

    pub fn pc_saft_da_res_ddens(&self, dens: f64, t: f64, x: &Array2<f64>) -> f64 {
        let h = 1e-4;
        let a_res_mais = self.pc_saft_a_res(dens + h, t, x);
        let a_res_menos = self.pc_saft_a_res(dens - h, t, x);
        let res = (a_res_mais - a_res_menos) / (2. * h);
        res
    }
    
    pub fn pc_saft_psat(&self, t: f64, guessp: f64) -> f64 {
        let x: Array2<f64> = Array2::from_elem((1,1), 1.0);

        let residuo = |var: f64| -> f64 {
            let psat = var*1e5;
            let (dens_l, dens_v) = self.pc_saft_dens(t, psat, &x, None);

            let phi_l = self.pc_saft_phi(dens_l, t, &x);
            let phi_v = self.pc_saft_phi(dens_v, t, &x);
            let res1 = 1.-phi_l[[0,0]]/phi_v[[0,0]];
            res1
        };

        let ans = find_interval(guessp.clone(), residuo);
        let psat = find_root(ans, residuo)*1e5;
        psat
    }

    pub fn pc_saft_massdens(&self, t: f64, p: f64, x: &Array2<f64>, phase: Option<&str>) -> Array2<f64> {
        let m2 = match &self.m2 {
            Some(m2_ref) => m2_ref.clone(),
            None => Array2::zeros((self.ncomp, 1)),
        };
        match phase {
            Some(k) => {
                let mut massdens: Array2<f64> = Array2::zeros((1, 1));
                let denslv = self.pc_saft_dens(t, p, x, Some(k));
                let dens = denslv.0 + denslv.1;
                for i in 0..self.ncomp {
                    massdens[[0, 0]] += x[[i, 0]] * m2[[i,0]] * dens / 1e+3;
                }
                massdens
            }
            None => {
                let mut massdens: Array2<f64> = Array2::zeros((2, 1));
                let (densl, densv) = self.pc_saft_dens(t, p, x, None);
                for i in 0..self.ncomp {
                    massdens[[0, 0]] += x[[i, 0]] * m2[[i,0]] * densl / 1e+3;
                    massdens[[1, 0]] += x[[i, 0]] * m2[[i,0]] * densv / 1e+3;
                }
                massdens
            }
        }
    }

    pub fn pc_saft_gamma_w(&self, dens: f64, t: f64, x: &Array2<f64>) -> Array2<f64> {
        let da_res_dx: Array2<f64> = self.pc_saft_da_res_dx(dens, t, x).unwrap();
        let mut gammaw: Array2<f64> = Array2::zeros((self.ncomp, 1));
        for i in 0..da_res_dx.len() {
            gammaw[[i, 0]] = f64::exp(da_res_dx[[i, 0]]);
        }
        gammaw
    }

    pub fn pc_saft_da_dt_num(&self, dens: f64, t: f64, x: &Array2<f64>) -> f64 {
        let step = 1e-5;
        
        let t_mais = t + step;
        let t_menos = t - step;

        let a_res_mais = self.pc_saft_a_res(dens, t_mais, x);
        let a_res_menos = self.pc_saft_a_res(dens, t_menos, x);
        let a_res_mais_mais = self.pc_saft_a_res(dens, t_mais + step, x);
        let a_res_menos_menos = self.pc_saft_a_res(dens, t_menos - step, x);

        let resp = (-a_res_mais_mais +8.*a_res_mais - 8.*a_res_menos+a_res_menos_menos)/(12.*step);
        resp
    }

    pub fn pc_saft_h_res_rt(&self, dens: f64, t: f64, x: &Array2<f64>) -> f64 {
        let da_dt = self.pc_saft_da_dt_num(dens, t, x);
        let z = self.pc_saft_z(dens, t, x);
        let h_res = -t*da_dt + (z-1.);
        h_res
    }

    pub fn pc_saft_s_res_r(&self, dens: f64, t: f64, x: &Array2<f64>) -> f64 {
        let a_res = self.pc_saft_a_res(dens, t, x);
        let da_dt = self.pc_saft_da_dt_num(dens, t, x);
        let z = self.pc_saft_z(dens, t, x);

        let s_res_r = -t*da_dt - a_res + f64::log(z, E); 
        s_res_r
    }

    pub fn pc_saft_g_res(&self, dens: f64, t: f64, x: &Array2<f64>) -> f64 {
        let a_res = self.pc_saft_a_res(dens, t, x);
        let z = self.pc_saft_z(dens, t, x);
        let g_res_rt = (a_res + (z-1.) - f64::log(z, E)) * KB * NAVO * t;
        g_res_rt
    }

    pub fn pc_saft_ddens_dt(&self, t: f64, p: f64, x: &Array2<f64>) -> Array2<f64> {
        let step = 1e-5;

        let t_mais = t + step;
        let t_menos = t - step;

        let dens_mais = self.pc_saft_massdens(t_mais, p, x, Some("vap"));
        let dens_menos = self.pc_saft_massdens(t_menos, p, x, Some("vap"));

        let ddens_dt: Array2<f64> = (dens_mais - dens_menos)/(2.*step);
        ddens_dt
    }

    pub fn pc_saft_alphap(&self, dens: f64, t: f64, x: &Array2<f64>) -> Array2<f64>{
        let ddens_dt = self.pc_saft_ddens_dt(t, dens, x);
        let dens = self.pc_saft_massdens(t, dens, x, Some("vap"));
        let alphap: Array2<f64> = -ddens_dt / dens;
        alphap
    }

    pub fn pc_saft_zassoc(&self, dens: f64, t: f64, x: &Array2<f64>) -> f64 {
        let mu_assoc_kt = self.pc_saft_mu_assoc_kt(dens, t, x).unwrap();
        let a_ass = self.pc_saft_a_ass(dens, t, x).unwrap();
        let mut zassoc: f64 = 0.;
        for i in 0..self.ncomp {
            zassoc += x[[i, 0]] * mu_assoc_kt[[i, 0]];
        }
        zassoc += a_ass;
        zassoc
    }

}
