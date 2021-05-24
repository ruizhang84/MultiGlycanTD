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
        public Dictionary<string, List<double>> Fragments { get; set; }
    }
}
