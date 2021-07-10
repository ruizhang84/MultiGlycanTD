using MultiGlycanTDLibrary.model;
using SpectrumProcess.algorithm;
using System;
using System.Collections.Generic;

namespace MultiGlycanTD
{
    public class SearchingParameters
    {
        // spectrum
        public double MS1Tolerance { get; set; } = 10;
        public double MSMSTolerance { get; set; } = 0.5;
        public ToleranceBy MS1ToleranceBy { get; set; } = ToleranceBy.PPM;
        public ToleranceBy MS2ToleranceBy { get; set; } = ToleranceBy.Dalton;

        // performance
        public int ThreadNums { get; set; } = 4;

        // searching
        public List<double> Ions { get; set; } = new List<double>()
            { MultiGlycanTDLibrary.util.mass.Spectrum.Proton };
        public double Similarity { get; set; } = 0.7;
        public double BinWidth { get; set; } = 1.0;

        // result
        public double FDR { get; set; } = 0.05;

        // file
        public List<string> MSMSFiles { get; set; } = new List<string>();
        public string DecoyFile { get; set; } = "";
        public GlycanJson Database { get; set; } = null;
        public string PeakFile { get; set; } = "";

        SearchingParameters() { }

        public void Update()
        {
            MS1Tolerance = ConfigureParameters.Access.MS1Tolerance;
            MSMSTolerance = ConfigureParameters.Access.MSMSTolerance;
            MS1ToleranceBy = ConfigureParameters.Access.MS1ToleranceBy;
            MS2ToleranceBy = ConfigureParameters.Access.MS2ToleranceBy;
            ThreadNums = ConfigureParameters.Access.ThreadNums;
            Similarity = ConfigureParameters.Access.Similarity;
            BinWidth = ConfigureParameters.Access.BinWidth;
            FDR = ConfigureParameters.Access.FDR;
            Ions = ConfigureParameters.Access.Ions;
        }

        protected static readonly Lazy<SearchingParameters>
            lazy = new Lazy<SearchingParameters>(() => new SearchingParameters());

        public static SearchingParameters Access { get { return lazy.Value; } }

    }

    public class ConfigureParameters
    {
        //spectrum
        public double MS1Tolerance { get; set; } = 10;
        public double MSMSTolerance { get; set; } = 0.5;
        public ToleranceBy MS1ToleranceBy { get; set; } = ToleranceBy.PPM;
        public ToleranceBy MS2ToleranceBy { get; set; } = ToleranceBy.Dalton;

        //Performance
        public int ThreadNums { get; set; } = 4;

        // searching
        public List<double> Ions { get; set; } = new List<double>()
            { MultiGlycanTDLibrary.util.mass.Spectrum.Proton };
        public double Similarity { get; set; } = 0.7;
        public double BinWidth { get; set; } = 1.0;

        //result
        public double FDR { get; set; } = 0.05;

        protected static readonly Lazy<ConfigureParameters>
            lazy = new Lazy<ConfigureParameters>(() => new ConfigureParameters());

        public static ConfigureParameters Access { get { return lazy.Value; } }

        private ConfigureParameters() { }

    }
}
