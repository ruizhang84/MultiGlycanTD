using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.model
{
    public class JsonEntry
    {
        public string ID { get; set; }
        public List<double> Fragments { get; set; }
    }

    public class GlycanJson
    {
       public CompdJson Compound { get; set; }
       public List<JsonEntry> Entries { get; set; }
    }
}
