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
    public class GlycanCompositionTest
    {
        [Test]
        public void Test()
        {
            IGlycan glycan = new NGlycanComplex();
            SortedDictionary<Monosaccharide, int> composition = new SortedDictionary<Monosaccharide, int>();
            composition[Monosaccharide.GlcNAc] = 2;
            composition[Monosaccharide.Man] = 3;
            glycan.SetComposition(composition);

            Glycan.To.SetPermethylation(true, true);
            Console.WriteLine( Glycan.To.Compute(glycan));

            GlycanBuilder builder = new GlycanBuilder();
            Compound compound = builder.BuildCompound(glycan);
            double mass = 0;
            foreach(var pair in compound.Composition)
            {
                ElementType ele = pair.Key;
                int count = pair.Value;

                mass += Brain.Run.ElementConv(ele).MonoMass * count;

            }
            Console.WriteLine(mass);

            Assert.AreEqual(Glycan.To.Compute(glycan), mass, 0.01);

        }

    }
}
