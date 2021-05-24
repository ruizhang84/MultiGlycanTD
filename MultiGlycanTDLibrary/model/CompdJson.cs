using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.model
{
    public class CompdJson
    {
        public Dictionary<string, List<string>> IDMap { get; set; }
        public Dictionary<string, List<double>> DistrMap { get; set; }
        public Dictionary<string, List<double>> MassMap { get; set; }
    }
}
