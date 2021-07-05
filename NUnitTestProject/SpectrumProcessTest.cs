using NUnit.Framework;
using SpectrumData;
using SpectrumData.Reader;
using SpectrumData.Spectrum;
using SpectrumProcess;
using SpectrumProcess.deisotoping;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NUnitTestProject
{
    public class SpectrumProcessTest
    {
        protected virtual List<IPeak> InsertPeaks(List<IPeak> origPeaks, 
            double precision=0.1)
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

        [Test]
        public void TestProcess()
        {
            // read spectrum
            string path = @"C:\Users\iruiz\Downloads\MSMS\HBS1_dextrinspkd_C18_10252018.raw";
            string output = @"C:\Users\iruiz\Downloads\MSMS\HBS1_dextrinspkd_C18_10252018_peaks.csv";
            ThermoRawSpectrumReader reader = new ThermoRawSpectrumReader();
            reader.Init(path);

            int start = reader.GetFirstScan();
            int end = reader.GetLastScan();
            IProcess process = new WeightedAveraging(new LocalNeighborPicking());
            //Averagine averagine = new Averagine(AveragineType.PermethylatedGlycan);
            //AveragineDeisotoping deisotoping = new AveragineDeisotoping(averagine,
            //    4, ToleranceBy.Dalton, 0.1);
            Dictionary<int, List<IPeak>> spectra = new Dictionary<int, List<IPeak>>();
            for (int i = start; i < end; i++)
            {
                if (reader.GetMSnOrder(i) == 2)
                {
                    ISpectrum ms2 = reader.GetSpectrum(i);
                    if (ms2.GetPeaks().Count <= 30)
                        continue;
                    ms2 = process.Process(ms2);
                    List<IPeak> ms2Peaks = ms2.GetPeaks();
                    spectra[i] = InsertPeaks(ms2Peaks);
                }
            }

            using (FileStream ostrm = new FileStream(output, FileMode.OpenOrCreate, FileAccess.Write))
            {
                using (StreamWriter writer = new StreamWriter(ostrm))
                {
                    writer.WriteLine("scan,mz,intensity");
                    foreach(int scan in spectra.OrderBy(s => s.Key).Select(s => s.Key))
                    {
                        List<IPeak> ms2Peaks = spectra[scan];
                        foreach (IPeak pk in ms2Peaks.OrderBy(p => p.GetMZ()))
                        {
                            writer.WriteLine(scan.ToString() + "," +
                                Math.Round(pk.GetMZ(), 3).ToString() + "," +
                                Math.Round(pk.GetIntensity(), 4).ToString());
                        }
                    }
                    
                }
            }
        }
    }
}