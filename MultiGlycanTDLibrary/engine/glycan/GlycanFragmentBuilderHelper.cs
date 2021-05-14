using MultiGlycanTDLibrary.model.glycan;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.glycan
{
    public class GlycanFragmentBuilderHelper
    {
     
        // assume glycan has GlcNAc(2)-man(3) full pentacore and at least 2 chains, at least 3 frgments
        public static int CountYCutComplex(IGlycan sub, IGlycan glycan, int limit)
        {
            int[] table = glycan.Table();
            int[] subTable = sub.Table();

            // at least 3 monsaccharides
            if (subTable[0] + subTable[1] + subTable[2] < 3)
                return -1;

            int diff = 0;
            // fucose core
            if (subTable[2] != table[2])
                diff++;
            // bisect
            if (subTable[3] != table[3])
                diff++;
            // not exceed limit
            if (diff > limit)
                return diff;

            // attenna not available
            if (subTable[1] == 0)
                return diff + 1;
            else if (subTable[4] == 0 && subTable[5] == 0)
                return diff + 2;
            else if (subTable[4] == 0 || subTable[5] == 0)
                return diff + 1;

            // branch
            for (int i = 0; i < 4; i++)
            {
                // branch
                if (table[i + 6] != subTable[i + 6] ||
                    table[i + 10] != subTable[i + 10] ||
                    table[i + 18] != subTable[i + 18] ||
                    table[i + 22] != subTable[i + 22])
                    diff++;
                // fucose only missing
                if (table[i + 6] == subTable[i + 6] &&
                    table[i + 14] != subTable[i + 14])
                    diff++;
                if (diff > limit)
                    break;
            }
            return diff;
        }

        public static int CountYCutHybrid(IGlycan sub, IGlycan glycan, int limit)
        {
            int[] table = glycan.Table();
            int[] subTable = sub.Table();

            // at least 3 monsaccharides
            if (subTable[0] + subTable[1] + subTable[2] < 3)
                return -1;

            int diff = 0;
            // fucose core
            if (subTable[2] != table[2])
                diff++;
            // bisect
            if (subTable[3] != table[3])
                diff++;
            // not exceed limit
            if (diff > limit)
                return diff;

            // attenna not available
            if (subTable[1] == 0)
                return diff + 1;
            else if (subTable[4] == 0 && subTable[5] == 0)
                return diff + 2;
            else if (subTable[4] == 0 || subTable[5] == 0)
                return diff + 1;

            // branch
            for (int i = 0; i < 2; i++)
            {
                if (table[i + 8] != subTable[i + 8] ||
                    table[i + 10] != subTable[i + 10] ||
                    table[i + 14] != subTable[i + 14] ||
                    table[i + 16] != subTable[i + 16])
                    diff++;
                if (table[i + 8] == subTable[i + 8] &&
                    table[i + 12] != subTable[i + 12])
                    diff++;
                if (table[i + 6] != subTable[i + 6])
                    diff++;
                if (diff > limit)
                    break;
            }
            return diff;
        }

        public static int CountYCutHighMannose(IGlycan sub, IGlycan glycan, int limit)
        {
            int[] table = glycan.Table();
            int[] subTable = sub.Table();

            // at least 3 monsaccharides
            if (subTable[0] + subTable[1] + subTable[2] < 3)
                return -1;

            int diff = 0;
            // fucose core
            if (subTable[2] != table[2])
                diff++;
            // not exceed limit
            if (diff > limit)
                return diff;

            // attenna not available
            if (subTable[1] == 0)
                return diff + 1;
            else if (subTable[3] == 0 && subTable[4] == 0)
                return diff + 2;
            else if (subTable[3] == 0 || subTable[4] == 0)
                return diff + 1;

            // branch
            for (int i = 0; i < 3; i++)
            {
                if (table[i + 5] != subTable[i + 5])
                    diff++;
                if (diff > limit)
                    break;
            }
            return diff;
        }

        public static int CountYCut(IGlycan sub, IGlycan glycan, int limit)
        {
            // branch
            switch (sub.Type())
            {
                case GlycanType.NGlycanComplex:
                    return CountYCutComplex(sub, glycan, limit);
                case GlycanType.NGlycanHybrid:
                    return CountYCutHybrid(sub, glycan, limit);
                case GlycanType.NHighMannose:
                    return CountYCutHighMannose(sub, glycan, limit);
            }
            return -1;
        }

        public static IGlycan ComplementaryFragment(IGlycan sub, IGlycan glycan)
        {
            IGlycan newGlycan = glycan.Clone();
            // compose
            SortedDictionary<Monosaccharide, int> compose = newGlycan.Composition();
            SortedDictionary<Monosaccharide, int> subCompose = sub.Composition();
            foreach(Monosaccharide sugar in compose.Keys)
            {
                compose[sugar] -= subCompose[sugar];
            }
            // core
            int[] table = newGlycan.Table();
            int[] subTable = sub.Table();
            for(int i = 0; i < table.Length; i++)
            {
                table[i] -= subTable[i];
            }
            // set glycan
            newGlycan.SetTable(table);
            newGlycan.SetComposition(compose);
            return newGlycan;
        }

    }
}
