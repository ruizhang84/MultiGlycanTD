using MultiGlycanClassLibrary.util.mass;
using MultiGlycanTDLibrary.engine.glycan;
using MultiGlycanTDLibrary.model.glycan;
using SpectrumProcess.brain;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NUnitTestProject
{
    public class GlycanMassCalculateUnitTest
    {
        [Test]
        public void GlycanTest()
        {
            // init
            Glycan.To.SetPermethylation(false, false);
            Glycan.To.Derivatization = Glycan.k2AB;
            GlycanBuilder glycanBuilder =
               new GlycanBuilder(7, 7, 5, 6, 0,
               true, false, false,
               10, false, false, Derivatization.k2AB);

            // composition
            IGlycan g = new NGlycanComplex();
            SortedDictionary<Monosaccharide, int> compos = new SortedDictionary<Monosaccharide, int>();
            compos[Monosaccharide.GlcNAc] = 4;
            compos[Monosaccharide.Man] = 3;
            compos[Monosaccharide.Gal] = 2;
            compos[Monosaccharide.Fuc] = 1;
            compos[Monosaccharide.NeuAc] = 1;
            g.SetComposition(compos);

            //Assert.AreEqual(2077.7349, Glycan.To.Compute(g), 0.1);


            // elements
            Dictionary<ElementType,int> compoundCompose = glycanBuilder.BuildCompound(g).Composition;
            double mass = 0;
            H h = new H();
            O o = new O();
            C c = new C();
            N n = new N();

            foreach(ElementType elem in compoundCompose.Keys)
            {
                switch(elem)
                {
                    case ElementType.C:
                        mass += compoundCompose[elem] * c.MonoMass;
                        break;
                    case ElementType.H:
                        mass += compoundCompose[elem] * h.MonoMass;
                        break;
                    case ElementType.O:
                        mass += compoundCompose[elem] * o.MonoMass;
                        break;
                    case ElementType.N:
                        mass += compoundCompose[elem] * n.MonoMass;
                        break;
                }    
            }

            Assert.AreEqual(Glycan.To.Compute(g), mass, 0.1);
        }
    }
}
