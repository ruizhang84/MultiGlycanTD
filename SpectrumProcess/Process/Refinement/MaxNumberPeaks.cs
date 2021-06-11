using SpectrumData;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SpectrumProcess
{
    public class MaxNumberPeaks: IProcess
    {
        IProcess process;
        int max;

        public MaxNumberPeaks(IProcess process, int max = 100)
        {
            this.process = process;
            this.max = max;
        }

        public List<IPeak> Process(List<IPeak> peaks)
        {
            List<IPeak> processed = process.Process(peaks);
            return processed
                .OrderByDescending(p => p.GetIntensity())
                .Take(max).ToList();
        }

        public ISpectrum Process(ISpectrum spectrum)
        {
            ISpectrum filtered = process.Process(spectrum);
            List<IPeak> processed = filtered.GetPeaks()
                .OrderByDescending(p => p.GetIntensity())
                .Take(max).ToList();
            filtered.SetPeaks(processed);
            return filtered;
        }
    }
}
