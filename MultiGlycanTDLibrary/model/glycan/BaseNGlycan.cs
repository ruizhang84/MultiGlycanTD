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
        protected HashSet<IGlycan> glycans = new HashSet<IGlycan>();
        protected SortedDictionary<Monosaccharide, int> composite
            = new SortedDictionary<Monosaccharide, int>();
        protected string id;
        protected object obj = new object();

        public abstract int[] Table();
        public abstract void SetTable(int[] table);

        public void Add(IGlycan glycan)
        {
            lock(obj)
            {
                glycans.Add(glycan);
            }
        }

        public List<IGlycan> Children()
        {
            return glycans.ToList();
        }

        public List<IGlycan> Fragments()
        {
            List<IGlycan> children = new List<IGlycan>();
            Stack<IGlycan> stack = new Stack<IGlycan>(glycans);
            HashSet<string> visited = new HashSet<string>();

            while(stack.Count > 0)
            {
                IGlycan node = stack.Pop();
                children.Add(node);

                foreach (IGlycan child in node.Children())
                {
                    string id = child.ID();
                    if(!visited.Contains(id))
                    {
                        visited.Add(id);
                        stack.Push(child);
                    }
                }
            }
            return children;
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

        public abstract GlycanType Type();
        public abstract IGlycan Clone();

        public virtual bool IsValid()
        {  
            return composite[Monosaccharide.GlcNAc] >= 2
                && composite[Monosaccharide.Man] >= 3;
        }

    }
}
