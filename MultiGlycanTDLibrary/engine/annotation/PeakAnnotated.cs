using MultiGlycanTDLibrary.engine.glycan;
using SpectrumData;
using System.Collections.Generic;

namespace MultiGlycanTDLibrary.engine.annotation
{
    public class GlycanAnnotated
    {
        public FragmentType Type { get; set; }
        public string Glycan { get; set; }
        public string Parent { get; set; }
    }
    public class PeakAnnotated
    {
        public IPeak Peak { get; set; }
        public string Glycan { get; set; }
        public List<GlycanAnnotated> Fragments { get; set; }
            = new List<GlycanAnnotated>();
    }
}
