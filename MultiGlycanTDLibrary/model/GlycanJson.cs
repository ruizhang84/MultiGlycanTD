using MultiGlycanTDLibrary.engine.glycan;
using System.Collections.Generic;

namespace MultiGlycanTDLibrary.model
{
    using GlycanFragments = Dictionary<FragmentType, List<string>>;

    public enum DerivationType
    {
        Permethylated, Native
    }

    public class GlycanJson
    {
        public DerivationType Derivation { get; set; }
        public CompdJson Compound { get; set; }
        // name(composition) -> id (structure)
        public Dictionary<string, List<string>> IDMap { get; set; }

        // fragments mass -> fragmenttype -> (intact/parent) glycan 
        public Dictionary<double, GlycanFragments> FragmentMap { get; set; }
        public ParameterJson Parameters { get; set; }
    }
}
