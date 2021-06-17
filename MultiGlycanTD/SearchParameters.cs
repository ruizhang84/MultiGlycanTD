using MultiGlycanTDLibrary.algorithm;
using MultiGlycanTDLibrary.model;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTD
{
    public class SearchingParameters
    {
        // spectrum
        public double MS1Tolerance { get; set; } = 10;
        public double MSMSTolerance { get; set; } = 0.1;
        public ToleranceBy MS1ToleranceBy { get; set; } = ToleranceBy.PPM;
        public ToleranceBy MS2ToleranceBy { get; set; } = ToleranceBy.Dalton;

        // performance
        public int ThreadNums { get; set; } = 4;

        // searching
        public List<double> Ions { get; set; } = new List<double>()
            { MultiGlycanTDLibrary.util.mass.Spectrum.Proton };

        // result
        public double Cutoff { get; set; } = 0.01;
        public double FDR { get; set; } = 0.01;

        // file
        public List<string> MSMSFiles { get; set; } = new List<string>();
        public GlycanJson Database { get; set; } = null;

        SearchingParameters()
        {
        }

        public void Update()
        {
            MS1Tolerance = ConfigureParameters.Access.MS1Tolerance;
            MSMSTolerance = ConfigureParameters.Access.MSMSTolerance;
            MS1ToleranceBy = ConfigureParameters.Access.MS1ToleranceBy;
            MS2ToleranceBy = ConfigureParameters.Access.MS2ToleranceBy;
            Cutoff = ConfigureParameters.Access.Cutoff;
            ThreadNums = ConfigureParameters.Access.ThreadNums;
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
        public double MSMSTolerance { get; set; } = 0.1;
        public ToleranceBy MS1ToleranceBy { get; set; } = ToleranceBy.PPM;
        public ToleranceBy MS2ToleranceBy { get; set; } = ToleranceBy.Dalton;

        //Performance
        public int ThreadNums { get; set; } = 4;

        // searching
        public List<double> Ions { get; set; } = new List<double>()
            { MultiGlycanTDLibrary.util.mass.Spectrum.Proton };

        //result
        public double Cutoff { get; set; } = 0.01;
        public double FDR { get; set; } = 0.01;

        protected static readonly Lazy<ConfigureParameters>
            lazy = new Lazy<ConfigureParameters>(() => new ConfigureParameters());

        public static ConfigureParameters Access { get { return lazy.Value; } }

        private ConfigureParameters() { }

    }
}
