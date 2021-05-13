using MultiGlycanTDLibrary.model.glycan;
using MultiGlycanTDLibrary.util.brain;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.glycan
{
    public class GlycanTheoryDistrbBuilder
    {
        bool permethylated = true;
        int order = 10;
        Dictionary<string, List<double>> mem;
        private readonly double neutron = 1.0;

        public GlycanTheoryDistrbBuilder(bool permethylated = true)
        {
            this.permethylated = permethylated;
            mem = new Dictionary<string, List<double>>();
        }

        public void SetPermethylated(bool permethylated)
        {
            this.permethylated = permethylated;
        }

        public IGlycan Build(IGlycan glycan)
        {
            Dictionary<Element, int> formulaComposition = new Dictionary<Element, int>();
            var compose = glycan.Composition();
            foreach (var sugar in compose.Keys)
            {
                Dictionary<Element, int> tempCompose = NMonosaccharideCreator.Get.SubCompositions(
                    sugar, permethylated);
                foreach (Element elm in tempCompose.Keys)
                {
                    if (!formulaComposition.ContainsKey(elm))
                    {
                        formulaComposition[elm] = 0;
                    }
                    formulaComposition[elm] += tempCompose[elm] * compose[sugar];
                }
            }
            Compound formula = new Compound(formulaComposition);
            glycan.SetFormula(formula);
            if (!mem.ContainsKey(formula.Name))
            {
                mem[formula.Name] = Brain.Run.Distribute(glycan.Formula(), order);
            }
            List<double> distrib = mem[formula.Name];

            glycan.SetDistrib(mem[formula.Name]);
            int extra = distrib.IndexOf(distrib.Max());
            glycan.SetHighestPeak(glycan.Mass() + neutron * extra);
            return glycan;
        }


        public Dictionary<string, IGlycan> Build(Dictionary<string, IGlycan> glycans_map)
        {
            foreach (string name in glycans_map.Keys)
            {
                IGlycan glycan = glycans_map[name];
                glycans_map[name] = Build(glycan);
            }
            return glycans_map;
        }

    }
}
