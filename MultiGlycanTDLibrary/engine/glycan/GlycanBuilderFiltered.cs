using MultiGlycanTDLibrary.model.glycan;
using MultiGlycanClassLibrary.util.mass;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MultiGlycanTDLibrary.util.brain;
using System.Collections.Concurrent;

namespace MultiGlycanTDLibrary.engine.glycan
{
    public class GlycanBuilderFiltered : GlycanBuilder
    {

        public GlycanBuilderFiltered(int hexNAc = 12, int hex = 12, int fuc = 5, int neuAc = 4, int neuGc = 0,
            bool complex = true, bool hybrid = false, bool highMannose = false, int order = 10,
            bool permethylated = true, bool reduced = true, int thread = 4): 
            base(hexNAc, hex, fuc, neuAc, neuGc, complex, 
                hybrid, highMannose, order, 
                permethylated, reduced, thread)
        {
        }

        public override Dictionary<string, IGlycan> GlycanMaps()
        { return glycans_map_; }

        protected override void BuildMaps()
        {
            foreach (var pair in glycans_map_)
            {
                IGlycan g = pair.Value;
                if (!g.IsValid())
                    continue;
                string name = g.Name();
                if (!compound_map_.ContainsKey(name))
                {
                    compound_map_[name] = BuildCompound(g);
                    glycan_compound_map_[name] = new List<IGlycan>();
                }
                glycan_compound_map_[name].Add(g);
            };

            Parallel.ForEach(compound_map_, pair =>
            {
                string name = pair.Key;
                Compound compound = pair.Value;
                distr_map_[name] = Brain.Run.Distribute(compound, Order);
                mass_map_[name] = Brain.Run.CenterMass(compound, Order);
            });
        }

    }
}
