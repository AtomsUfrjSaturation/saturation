using System;
using HYSYS;
using System.Runtime.InteropServices;
using System.IO;

namespace Atoms
{
    public delegate void WaterFraction(double data);

    [ComVisible(true)]
    [ProgId("Atoms.Saturation")]
    [ClassInterface(ClassInterfaceType.AutoDispatch)]
    public class Saturation
    {
        private void Logger(string message)
        {
            string desktopPath = Environment.GetFolderPath(Environment.SpecialFolder.Desktop);
            string extensionFolder = Path.Combine(desktopPath, "ext");  // change this depending on where the DLL is
            string filePath = Path.Combine(extensionFolder, "csharp_log.txt");

            string timestamp = DateTime.UtcNow.ToString("yyyy-MM-dd HH:mm:ss.fff");
            string content = $"{timestamp}: {message}\n";

            File.AppendAllText(filePath, content);
        }


        [DllImport("atoms_saturation_kernel.dll")]
        public static extern void send_info_to_hysys(WaterFraction callback);

        [DllImport("atoms_saturation_kernel.dll")]
        public static extern void receive_info_from_hysys(IntPtr value, IntPtr label);

        private void OnWaterFractionReceived(double value)
        {
            string message = $"Received water fraction from the modelling libray: {value}\n";
            Logger(message);
            waterFraction.SetValue(value);
        }

        public void GetWaterFraction()
        {
            send_info_to_hysys(OnWaterFractionReceived);
        }

        private void SendToHysys(string temperature, string pressure)
        {
            IntPtr temperaturePtr = Marshal.StringToHGlobalAnsi(temperature);
            IntPtr pressurePtr = Marshal.StringToHGlobalAnsi(pressure);
            try
            {
                receive_info_from_hysys(temperaturePtr, pressurePtr);
            }
            finally
            {
                Marshal.FreeHGlobal(temperaturePtr);
                Marshal.FreeHGlobal(pressurePtr);
            }
        }

        #region: Class Variables declaration
        private HYSYS.ExtnUnitOperationContainer extnContainer;
        private HYSYS.ProcessStream Feed;
        private HYSYS.ProcessStream Product;
        private HYSYS.InternalRealVariable eosModel;
        private HYSYS.InternalRealVariable assocModel;
        private HYSYS.InternalRealVariable waterFraction;
        private HYSYS.InternalRealVariable crossH2S;
        private HYSYS.InternalRealVariable crossCO2;
        const int emptyVal = (int)EmptyValue_enum.HEmpty;
        #endregion

        public int Initialize(HYSYS.ExtnUnitOperationContainer Container, bool IsRecalling)
        {
            extnContainer = Container;
            Feed = extnContainer.FindVariable("FeedStream").Variable.Object;
            Product = extnContainer.FindVariable("ProductStream").Variable.Object;
            eosModel = extnContainer.FindVariable("Eos_Model").Variable;
            assocModel = extnContainer.FindVariable("Assoc_Model").Variable;
            waterFraction = extnContainer.FindVariable("Water_Fraction").Variable;
            crossH2S = extnContainer.FindVariable("Cross_H2S_sol").Variable;
            crossCO2 = extnContainer.FindVariable("Cross_CO2_sol").Variable;

            return (int)HYSYS.CurrentExtensionVersion_enum.extnCurrentVersion;
        }

        public void Execute(bool isForgetpass)
        {
            try
            {
                if (isForgetpass) return;

                Feed = extnContainer.FindVariable("FeedStream").Variable.Object;
                Product = extnContainer.FindVariable("ProductStream").Variable.Object;
                eosModel = extnContainer.FindVariable("Eos_Model").Variable;
                assocModel = extnContainer.FindVariable("Assoc_Model").Variable;
                waterFraction = extnContainer.FindVariable("Water_Fraction").Variable;
                crossH2S = extnContainer.FindVariable("Cross_H2S_sol").Variable;
                crossCO2 = extnContainer.FindVariable("Cross_CO2_sol").Variable;


                if (Feed != null &&
                    Feed.Temperature.IsKnown &&
                    Feed.Pressure.IsKnown)
                {
                    double temp = Feed.Temperature.GetValue();
                    double pres = Feed.Pressure.GetValue();

                    if (temp == 0 || pres == 0) return;

                    SendToHysys(temp.ToString(), pres.ToString());

                    try
                    {
                        GetWaterFraction();
                    }
                    catch (Exception ex)
                    {
                        Logger(ex.ToString());
                    }
                }
            }
            catch (Exception ex)
            {
                Logger(ex.ToString());
            }
        }
    }
}
