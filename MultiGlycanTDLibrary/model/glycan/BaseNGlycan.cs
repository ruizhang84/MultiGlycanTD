using MultiGlycanTDLibrary.util.brain;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.model.glycan
{
    public abstract class BaseNGlycan : IGlycan
    {
        protected Dictionary<string, IGlycan> glycans = new Dictionary<string, IGlycan>();
        protected SortedDictionary<Monosaccharide, int> composite
            = new SortedDictionary<Monosaccharide, int>();
        protected string id;
        protected double mass;

        // get distribution of isotope cluster
        protected Compound formula;
        protected List<double> distrib;
        protected double peak = 0;

        public abstract int[] Table();
        public abstract void SetTable(int[] table);

        public void Add(IGlycan glycan)
        {
            glycans[glycan.ID()] = glycan;
            foreach (IGlycan g in glycan.Children())
            {
                glycans[g.ID()] = g;
            }
        }

        public List<IGlycan> Children()
        {
            return glycans.Values.ToList();
        }

        public SortedDictionary<Monosaccharide, int> Composition()
        {
            return composite;
        }

        public abstract List<IGlycan> Grow(Monosaccharide monosaccharide);

        public string ID()
        {
            return string.Join(" ", Table()); ;
        }

        public double Mass()
        {
            return mass;
        }

        public string Name()
        {
            string name = "";
            foreach (Monosaccharide sugar in composite.Keys)
            {
                switch (sugar)
                {
                    case Monosaccharide.GlcNAc:
                        name += "GlcNAc-" + composite[sugar] + " ";
                        break;
                    case Monosaccharide.Man:
                        name += "Man-" + composite[sugar] + " ";
                        break;
                    case Monosaccharide.Gal:
                        name += "Gal-" + composite[sugar] + " ";
                        break;
                    case Monosaccharide.Fuc:
                        name += "Fuc-" + composite[sugar] + " ";
                        break;
                    case Monosaccharide.NeuAc:
                        name += "NeuAc-" + composite[sugar] + " ";
                        break;
                    case Monosaccharide.NeuGc:
                        name += "NeuGc-" + composite[sugar] + " ";
                        break;
                    default:
                        break;
                }
            }
            return name;
        }

        public void SetComposition(SortedDictionary<Monosaccharide, int> composite)
        {
            this.composite = new SortedDictionary<Monosaccharide, int>();
            foreach (var key in composite.Keys)
            {
                this.composite[key] = composite[key];
            }
        }

        public void SetMass(double mass)
        {
            this.mass = mass;
        }

        public abstract GlycanType Type();
        public abstract IGlycan Clone();

        public virtual bool IsValid()
        {  
            return composite[Monosaccharide.GlcNAc] >= 2
                && composite[Monosaccharide.Man] >= 3;
        }

        public Compound Formula()
        {
            return formula;
        }

        public void SetFormula(Compound formula)
        {
            this.formula = formula;
        }

        public List<double> GetDistrib()
        {
            return distrib;
        }

        public void SetDistrib(List<double> distrib)
        {
            this.distrib = distrib;
        }

        public double HighestPeak()
        {
            return peak;
        }
        public void SetHighestPeak(double peak)
        {
            this.peak = peak;
        }

    }
}
