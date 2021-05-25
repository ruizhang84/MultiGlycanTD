using MultiGlycanTDLibrary.engine.analysis;
using SpectrumData;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTD
{
    public class ProgressingEventArgs : EventArgs
    {
        public int Total { get; set; }
    }

    public class Counter
    {
        public event EventHandler<ProgressingEventArgs> progressChange;

        protected virtual void OnProgressChanged(ProgressingEventArgs e)
        {
            EventHandler<ProgressingEventArgs> handler = progressChange;
            handler?.Invoke(this, e);
        }

        public void Add(int total)
        {
            ProgressingEventArgs e = new ProgressingEventArgs
            {
                Total = total
            };
            OnProgressChanged(e);
        }
    }

    public class SearchTask
    {
        public SearchTask(ISpectrum spectrum, List<IPeak> peaks, double mz, int charge)
        {
            Spectrum = spectrum;
            MSPeaks = peaks;
            PrecursorMZ = mz;
            Charge = charge;
        }
        public ISpectrum Spectrum { get; set; }
        public List<IPeak> MSPeaks { get; set; }
        public double PrecursorMZ { get; set; }
        public int Charge { get; set; }

    }

    public class MultiThreadingSearchHelper
    {

        public static List<IPeak> FilterPeaks(List<IPeak> peaks, double target, double range)
        {
            if (peaks.Count == 0)
            {
                return peaks;
            }

            int start = 0;
            int end = peaks.Count - 1;
            int middle = 0;
            if (peaks[start].GetMZ() > target - range)
            {
                middle = start;
            }
            else
            {
                while (start + 1 < end)
                {
                    middle = (end - start) / 2 + start;
                    double mz = peaks[middle].GetMZ() + range;
                    if (mz == target)
                    {
                        break;
                    }
                    else if (mz < target)
                    {
                        start = middle;
                    }
                    else
                    {
                        end = middle - 1;
                    }
                }
            }

            List<IPeak> res = new List<IPeak>();
            while (middle < peaks.Count)
            {
                if (peaks[middle].GetMZ() > target + range)
                    break;
                res.Add(peaks[middle++]);
            }
            return res;
        }
        public static void ReportResults(
            string path, 
            List<SearchResult> results)
        {
            using (FileStream ostrm = new FileStream(path, FileMode.OpenOrCreate, FileAccess.Write))
            {
                using (StreamWriter writer = new StreamWriter(ostrm))
                {
                    writer.WriteLine("scan,retention,glycan,precursor_mz,score,fit");
                    foreach (SearchResult r in results.OrderBy(p => p.Scan()))
                    {
                        string output = r.Scan().ToString() + ","
                            + r.Retention().ToString() + ","
                            + r.Glycan() + ","
                            + r.MZ().ToString() + ","
                            + r.Score().ToString() + ","
                            + r.Fit().ToString();
                        writer.WriteLine(output);
                    }
                    writer.Flush();
                }
            }
        }

    }
}
