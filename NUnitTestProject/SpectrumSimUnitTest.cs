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
    public class SpectrumSimUnitTest
    {
        [Test]
        public void CosTest()
        {
            // read spectrum
            string path = @"C:\Users\iruiz\Downloads\MSMS\134144_31_C18_120min_60oC_50cm.raw";
            ThermoRawSpectrumReader reader = new ThermoRawSpectrumReader();
            reader.Init(path);

            ISpectrum A = reader.GetSpectrum(18512);
            ISpectrum B = reader.GetSpectrum(18746);

            Console.WriteLine(GlycanScorerHelper.CosineSim(A.GetPeaks(), B.GetPeaks(), 1.0));
                   
            Assert.Pass();
        }

     
    }
}
