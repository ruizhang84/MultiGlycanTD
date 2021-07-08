using SpectrumData;
using SpectrumData.Spectrum;
using System.Collections.Generic;
using System.Linq;

namespace SpectrumProcess
{
    public class LocalNeighborPicking : IProcess
    {
        protected double precision = 0.1;

        public LocalNeighborPicking(double precision = 0.1)
        {
            this.precision = precision;
        }

        protected virtual List<IPeak> InsertPeaks(List<IPeak> origPeaks)
        {
            List<IPeak> peaks = new List<IPeak>();
            double last = origPeaks.First().GetMZ();
            peaks.Add(new GeneralPeak(last - precision, 0));
            foreach (IPeak peak in origPeaks)
            {
                if (peak.GetMZ() - last > precision)
                {
                    peaks.Add(new GeneralPeak(last + precision / 2, 0));
                    peaks.Add(new GeneralPeak(peak.GetMZ() - precision / 2, 0));
                }
                peaks.Add(peak);
                last = peak.GetMZ();
            }
            peaks.Add(new GeneralPeak(last + precision, 0));
            return peaks;
        }

        public List<IPeak> Process(List<IPeak> peaks)
        {
            // insert pseudo peaks for large gap
            peaks = InsertPeaks(peaks);
            List<IPeak> processedPeaks = new List<IPeak>();

            int index = 1;
            int end = peaks.Count - 1;
            int head = index + 1;
            while (index < end)
            {
                if (peaks[index - 1].GetIntensity() < peaks[index].GetIntensity())
                {
                    head = index + 1;
                }

                while (head < end
                    && peaks[head].GetIntensity() == peaks[index].GetIntensity())
                {
                    head++;
                }

                if (peaks[head].GetIntensity() < peaks[index].GetIntensity())
                {
                    processedPeaks.Add(peaks[index]);
                    index = head;
                }
                index++;

            }
            return processedPeaks;
        }


        public ISpectrum Process(ISpectrum spectrum)
        {
            if (spectrum.GetPeaks().Count == 0)
                return spectrum;
            // insert pseudo peaks for large gap
            List<IPeak> peaks = InsertPeaks(spectrum.GetPeaks());
            List<IPeak> processedPeaks = Process(peaks);

            ISpectrum newSpectrum = spectrum.Clone();
            newSpectrum.SetPeaks(processedPeaks);
            return newSpectrum;
        }
    }
}
