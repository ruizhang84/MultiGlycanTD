using System.Collections.Generic;

namespace MultiGlycanTDLibrary.model
{
    public class CompdJson
    {
        public Dictionary<string, List<double>> DistrMap { get; set; }
        public Dictionary<string, List<double>> MassMap { get; set; }
    }
}
