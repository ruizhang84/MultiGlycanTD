using SpectrumData;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SpectrumProcess
{
    public class PeakNormalization : IProcess
    {
        public IProcess peakPicking;
        public PeakNormalization(IProcess peakPicking)
        {
            this.peakPicking = peakPicking;
        }

        public List<IPeak> Process(List<IPeak> peaks)
        {
            List<IPeak> processed = peakPicking.Process(peaks);
            double sum = processed.Select(p => Math.Sqrt(p.GetIntensity())).Sum();
            foreach (IPeak pk in processed)
            {
                pk.SetIntensity(Math.Sqrt(pk.GetIntensity()) / sum);
            }
            return processed;
        }

        public ISpectrum Process(ISpectrum spectrum)
        {
            ISpectrum processed = peakPicking.Process(spectrum);
            List<IPeak> peaks = Process(processed.GetPeaks());
            processed.SetPeaks(peaks);
            return processed;
        }
    }
}
