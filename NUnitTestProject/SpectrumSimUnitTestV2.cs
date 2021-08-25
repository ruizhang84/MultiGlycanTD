using SpectrumProcess.algorithm;
using MultiGlycanTDLibrary.engine.score;
using MultiGlycanTDLibrary.engine.search;
using MultiGlycanTDLibrary.model;
using NUnit.Framework;
using SpectrumData;
using SpectrumData.Reader;
using SpectrumProcess;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using System.Threading.Tasks;

namespace NUnitTestProject
{
    public class SpectrumSimUnitTestV2
    {
        [Test]
        public void CosTestV2()
        {
            // read spectrum
            string path = @"C:\Users\iruiz\Downloads\MultiGlycanTD\data\glycan_standard3.mgf";
            string output = @"C:\Users\iruiz\Downloads\MultiGlycanTD\data\glycan_standard3.csv";
            MGFSpectrumReader reader = new MGFSpectrumReader();
            reader.Init(path);

            using (FileStream ostrm = new FileStream(output, FileMode.OpenOrCreate, FileAccess.Write))
            {
                using (StreamWriter writer = new StreamWriter(ostrm))
                {
                    writer.WriteLine("scan1,scan2,cosine");
                    foreach (int i in reader.GetSpectrum().Keys)
                    {
                        foreach (int j in reader.GetSpectrum().Keys)
                        {
                            if (i == j) continue;
                            ISpectrum A = reader.GetSpectrum(i);
                            ISpectrum B = reader.GetSpectrum(j);
                            double cosine = GlycanScorerHelper.CosineSim(A.GetPeaks(), B.GetPeaks(), 1.0);
                            writer.WriteLine(i.ToString() + "," + j.ToString() + "," + cosine.ToString());
                        }
                    }

                }
            }

            Assert.Pass();
        }

     
    }
}
