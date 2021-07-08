using System.Collections.Generic;

namespace SpectrumData.Spectrum
{
    public class MS2Spectrum : ISpectrum
    {
        protected int scanNum;
        protected double retention;
        protected List<IPeak> peaks;
        protected double precursorMZ;
        protected int precursorCharge;
        protected TypeOfMSActivation activation;

        public MS2Spectrum(int scanNum, double retention,
            double precursorMZ, int precursorCharge)
        {
            this.scanNum = scanNum;
            this.retention = retention;
            this.precursorMZ = precursorMZ;
            this.precursorCharge = precursorCharge;
            peaks = new List<IPeak>();
        }


        public MS2Spectrum(ISpectrum spectrum,
            double mz, int charge)
        {
            scanNum = spectrum.GetScanNum();
            retention = spectrum.GetRetention();
            precursorMZ = mz;
            precursorCharge = charge;
            peaks = spectrum.GetPeaks();
        }

        public TypeOfMSActivation Activation()
        {
            return activation;
        }

        public void SetActivation(TypeOfMSActivation type)
        {
            activation = type;
        }

        public double PrecursorMZ()
        {
            return precursorMZ;
        }

        public int PrecursorCharge()
        {
            return precursorCharge;
        }

        public void Add(IPeak peak)
        {
            peaks.Add(peak);
        }
        public void Clear()
        {
            peaks.Clear();
        }

        public ISpectrum Clone()
        {
            ISpectrum spec = new MS2Spectrum(scanNum, retention,
                precursorMZ, precursorCharge);
            spec.SetPeaks(peaks);
            return spec;
        }

        public List<IPeak> GetPeaks()
        {
            return peaks;
        }

        public double GetRetention()
        {
            return retention;
        }

        public int GetScanNum()
        {
            return scanNum;
        }

        public void SetPeaks(List<IPeak> peaks)
        {
            this.peaks = peaks;
        }

        public void set_scan(int scan)
        {
            scanNum = scan;
        }

        public void set_retention(double retention)
        {
            this.retention = retention;
        }
    }
}
