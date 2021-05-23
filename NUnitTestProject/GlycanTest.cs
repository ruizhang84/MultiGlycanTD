using MultiGlycanTDLibrary.engine.glycan;
using MultiGlycanTDLibrary.model.glycan;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NUnitTestProject
{
    public class GlycanTest
    {
        [Test]
        public void Test()
        {
            int[] table = new int[26]
            {
                //2, 1, 0, 0, 1, 1, 2, 2, 2, 0, 2, 2, 2, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0,
                2, 1, 0, 0, 1, 1, 2, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            };

            IGlycan glycan = new NGlycanComplex();
            glycan.SetTable(table);
            List<IGlycan> gs = glycan.Grow(Monosaccharide.Fuc);
            //2 1 0 0 1 1 2 2 2 0 2 2 2 0 1 1 1 0 1 1 0 0 0 0 0 0
            // 2 1 0 0 1 1 2 2 2 0 2 2 2 0 1 1 1 0 1 0 1 0 0 0 0 0
            foreach (var g in gs)
            {
                Console.WriteLine(g.ID());
            }

            //GlycanBuilder glycanBuilder =
            //   new GlycanBuilder(12, 12, 5, 4, 0, true, false, false);
            //glycanBuilder.SpeedUp = true;
            //glycanBuilder.Build();

            ////var map = glycanBuilder.GlycanMaps();
            ////Assert.IsTrue(map.ContainsKey("2 1 0 0 1 1 2 2 2 0 2 2 2 0 1 1 1 0 1 1 0 0 0 0 0 0"));
            //Assert.IsTrue(glycanBuilder.SatisfyCriteria(glycan));

            Assert.Pass();

        }
    }
}
