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
        public void CosTest()
        {
            // read spectrum
            string path = @"C:\Users\iruiz\Downloads\MSMS\HBS1_dextrinspkd_C18_10252018.raw";
            ThermoRawSpectrumReader reader = new ThermoRawSpectrumReader();
            reader.Init(path);

            ISpectrum A = reader.GetSpectrum(930);
            ISpectrum B = reader.GetSpectrum(974);

            Console.WriteLine(GlycanScorerHelper.CosineSim(A.GetPeaks(), B.GetPeaks(), 1.0));
                   
            Assert.Pass();
        }

     
    }
}
