using System.Collections.Generic;
using System.Linq;

namespace SpectrumData.Spectrum
{
    public class GeneralSpectrum : ISpectrum
    {
        protected int scanNum;
        protected double retention;
        protected List<IPeak> peaks;

        public GeneralSpectrum(int scanNum, double retention)
        {
            this.scanNum = scanNum;
            this.retention = retention;
            peaks = new List<IPeak>();
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
            ISpectrum spec = new GeneralSpectrum(scanNum, retention);
            spec.SetPeaks(peaks.Count > 0 ?
                peaks.Select(p => p.Clone()).ToList() : new List<IPeak>());
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

        public void set_retention(double retention)
        {
            this.retention = retention;
        }

        public void set_scan(int scan)
        {
            scanNum = scan;
        }
    }
}
