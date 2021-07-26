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
            string path = @"C:\Users\iruiz\Downloads\MSMS\data_analysis\glycan_sample\Fetuin_1ug_01_C18_50cm_091620.mgf";
            string path2 = @"C:\Users\iruiz\Downloads\MSMS\data_analysis\glycan_sample\Fetuin_1ug_02_C18_50cm_091620.mgf";
            ISpectrumReader reader = new MGFSpectrumReader();
            reader.Init(path);
            ISpectrum A = reader.GetSpectrum(8160);

            reader.Init(path2);
            ISpectrum B = reader.GetSpectrum(7988);

            Console.WriteLine(GlycanScorerHelper.CosineSim(A.GetPeaks(), B.GetPeaks(), 1.0));
                   
            Assert.Pass();
        }

     
    }
}
