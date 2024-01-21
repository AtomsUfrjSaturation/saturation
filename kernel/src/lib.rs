mod pc_saft;

use std::fs::OpenOptions;
use std::io::Write;
use pc_saft::PcSaftEos;
use ndarray::Array2;
use std::sync::Mutex;
use std::os::raw::c_double;
use chrono::Utc; 

fn logger(message: &str) -> std::io::Result<()> {
    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .write(true)
        .open("rust_log.txt")?;

    let timestamp = Utc::now().format("%Y-%m-%d %H:%M:%S%.3f");
    let formatted_string = format!("{}: {}\n", timestamp, message);
    file.write_all(formatted_string.as_bytes())?;

    Ok(())
}

struct HysysProperties<'a> {
    temperature: &'a str, 
    pressure: &'a str, 
    //components: &'a str
}

static hysysProperties: Mutex<HysysProperties> = Mutex::new(
    HysysProperties {
        temperature: "",
        pressure: "",
        //components: ""
    }
);


#[no_mangle]
pub extern "C" fn receive_info_from_hysys(
    temperatureValue: *const libc::c_char,
    pressureValue: *const libc::c_char) -> () {
    
    let temperature_c_str = unsafe {
        std::ffi::CStr::from_ptr(temperatureValue)
    };
    let temperature_string = temperature_c_str.to_str().unwrap();

    let pressure_c_str = unsafe {
        std::ffi::CStr::from_ptr(pressureValue)
    };
    let pressure_string = pressure_c_str.to_str().unwrap();

    {
        let mut properties = hysysProperties.lock().unwrap();

        properties.temperature = temperature_string;
        properties.pressure = pressure_string;
    }
    perform_calculations();
}


fn perform_calculations() -> f64  {
    let properties = hysysProperties.lock().unwrap();
    let temperature = properties.temperature.parse::<f64>().unwrap_or_default();
    let pressure = properties.pressure.parse::<f64>().unwrap_or_default();

    let log_message= format!("Temperature: {}, Pressure: {}\n", temperature, pressure);
    let _ = logger(&log_message);

    // The block below contains artificial data for creating an 
    // instance of the PC-SAFT EoS struct.
    // TODO: harcdode the values of components properties
    // or get them from Hysys?
    let mut m1: Array2<f64> = Array2::zeros((1,1));
    m1[[0,0]] = 3.01016;

    let mut sigma: Array2<f64> = Array2::zeros((1,1));
    sigma[[0,0]] = 2.969;

    let mut epsilon_k: Array2<f64> = Array2::zeros((1, 1));
    epsilon_k[[0,0]] = 236.539;

    let mut kab_k: Array2<f64> = Array2::zeros((1, 1));
    kab_k[[0,0]] = 0.04;

    let mut kbi = Array2::zeros((1,1));
    kbi[[0,0]] = 0.;

    let deltasimplified = true;

    let eos = PcSaftEos::new(m1, sigma, epsilon_k, None, kbi, None, None, None, deltasimplified);
    let x = Array2::zeros((1,1));

    let data = eos.pc_saft_massdens(temperature, pressure, &x, None);
    // In the meantime, send a random hardcoded value for water fraction 
    // Because of the FFI between C# and Rust, this value will show up in 
    // the extension's view in Hysys
    0.3953105
}

#[no_mangle]
pub extern "C" fn send_info_to_hysys(callback: extern "C" fn(c_double)) {
    let waterFrac = perform_calculations();
    callback(waterFrac);
}
