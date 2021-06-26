using MultiGlycanTDLibrary.engine.glycan;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.model
{

    public class GlycanJson
    {
        public CompdJson Compound { get; set; }
        // name(composition) -> id (structure)
        public Dictionary<string, List<string>> IDMap { get; set; }

        // fragments mass -> fragmenttype -> (intact/parent) glycan 
        public Dictionary<double, Dictionary<FragmentTypes, List<string>>> FragmentMap { get; set; }
    }
}
