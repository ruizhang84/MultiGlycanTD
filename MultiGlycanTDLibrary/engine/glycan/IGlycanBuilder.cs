using MultiGlycanTDLibrary.model.glycan;
using System.Collections.Generic;

namespace MultiGlycanTDLibrary.engine.glycan
{
    public interface IGlycanBuilder
    {
        Dictionary<string, IGlycan> GlycanMaps();
        Dictionary<string, List<IGlycan>> GlycanCompositionMaps();
        Dictionary<string, List<double>> GlycanDistribMaps();
        Dictionary<string, List<double>> GlycanMassMaps();
        void Build();
    }
}
