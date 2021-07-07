using MultiGlycanTDLibrary.engine.glycan;
using SpectrumData;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.search
{
    public class PeakMatch
    {
        // matched peak
        public IPeak Peak { get; set; }
        // difference between theortical and observed peaks m/z
        public double TheoreticMZ { get; set; } = int.MaxValue;
        // potential matched to other glycans
        public int Potentials { get; set; } = int.MaxValue;
        // fragment ion types
        public HashSet<FragmentType> IonTypes { get; set; }
            = new HashSet<FragmentType>();
    }

    public class SearchResult
    {
        public int Scan { get; set; }
        public double Retention { get; set; }
        public double PrecursorMZ { get; set; }
        public double Ion { get; set; }
        public int Charge { get; set; }
        public string Glycan { get; set; }
        public string Composition { get; set; }
        public Dictionary<int, PeakMatch> Matches
            = new Dictionary<int, PeakMatch>();
        public double Score { get; set; }
        public double Fit { get; set; }
        public double Coverage { get; set; }
    }
}
