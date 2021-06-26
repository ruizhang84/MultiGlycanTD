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
        public double Diff { get; set; } = int.MaxValue;
        // potential matched to other glycans
        public int Potentials { get; set; } = int.MaxValue;
        // fragment ion types
        public HashSet<FragmentTypes> IonTypes { get; set; }
            = new HashSet<FragmentTypes>();
    }

    public class SearchResult
    {
        public int Scan { get; set; }
        public string Glycan { get; set; }
        public string Composite { get; set; }
        public Dictionary<int, PeakMatch> Matches
            = new Dictionary<int, PeakMatch>();
    }
}
