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

            int diff = 0;

            // partial penta core
            if (subTable[0] + subTable[1] + subTable[2] < 3)
            {
                if (subTable[0] == 0)
                    return -1;
                // fucose core
                if (subTable[2] != table[2])
                    diff++;
                diff++;
                return diff;
            }
            
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
                diff++;

            // branch sorted only for 1 and 2 cleavages
            if (glycan.Sorted())
            {
                // find equal chain
                int subIndex = 0;
                int index = 0;
                List<int> equalIndex = new List<int>();
                List<int> equalSubIndex = new List<int>();
                while (subIndex < 4 && index < 4)
                {
                    int length = table[index + 6] + table[index + 10] + table[index + 18];
                    int subLength = table[subIndex + 6] + table[subIndex + 10] + table[subIndex + 18];

                    if (length == subLength && table[index + 14] == subTable[subIndex + 14])
                    {
                        equalIndex.Add(index);
                        equalSubIndex.Add(subIndex);
                        subIndex++;
                        index++;
                    }
                    else if (length == subLength && table[index + 14] > subTable[subIndex + 14])
                    {
                        index++;
                    }
                    else if (length == subLength && table[index + 14] < subTable[subIndex + 14])
                    {
                        subIndex++;
                    }
                    else if (length > subLength)
                    {
                        index++;
                    }
                    else
                    {
                        subIndex++;
                    }
                }
                // enumerate and compute diff
                int startIndex = 0;
                int startSubIndex = 0;
                List<Tuple<int, int>> pairs = new List<Tuple<int, int>>();
                while (startIndex < 4 && startSubIndex < 4)
                {
                    if (equalIndex.Contains(startIndex))
                    {
                        startIndex++;
                        continue;
                    }
                    if (equalSubIndex.Contains(startSubIndex))
                    {
                        startSubIndex++;
                        continue;
                    }
                    pairs.Add(new Tuple<int, int>(startIndex, startSubIndex));
                    startIndex++;
                    startSubIndex++;
                }

                foreach(Tuple<int,int> indexPair in pairs)
                {
                    // branch
                    int i = indexPair.Item1;
                    int j = indexPair.Item2;
                    if(subTable[5] == 0 && j >= 2)
                    {
                        break;
                    }
                    if (table[i + 6] != subTable[j + 6] ||
                        table[i + 10] != subTable[j + 10] ||
                        table[i + 18] != subTable[j + 18] ||
                        table[i + 22] != subTable[j + 22])
                        diff++;
                    // fucose only missing
                    if (table[i + 6] == subTable[j + 6] &&
                        table[i + 14] != subTable[j + 14])
                        diff++;
                }
            }
            else
            {
                for (int i = 0; i < 2; i++)
                {
                    if (subTable[i + 4] > 0)
                    {
                        for (int j = 0; j < 2; j++)
                        {
                            // branch
                            if (table[i * 2 + j + 6] != subTable[i * 2 + j + 6] ||
                                table[i * 2 + j + 10] != subTable[i * 2 + j + 10] ||
                                table[i * 2 + j + 18] != subTable[i * 2 + j + 18] ||
                                table[i * 2 + j + 22] != subTable[i * 2 + j + 22])
                                diff++;
                            // fucose only missing
                            if (table[i * 2 + j + 6] == subTable[i * 2 + j + 6] &&
                                table[i * 2 + j + 14] != subTable[i * 2 + j + 14])
                                diff++;
                        }
                        if (diff > limit)
                            break;
                    }
                }
            }

            return diff;
        }

        public static int CountYCutHybrid(IGlycan sub, IGlycan glycan, int limit)
        {
            int[] table = glycan.Table();
            int[] subTable = sub.Table();

            int diff = 0;
            // partial fivecore
            if (subTable[0] + subTable[1] + subTable[2] < 3)
            {
                if (subTable[0] == 0)
                    return -1;
                // fucose core
                if (subTable[2] != table[2])
                    diff++;
                diff++;
                return diff;
            }

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
                diff++;

            // branch
            if (subTable[4] > 0)
            {
                for (int i = 0; i < 2; i++)
                {
                    if (table[i + 6] != subTable[i + 6])
                        diff++;
                    if (diff > limit)
                        break;
                }
            }
            if (subTable[5] > 0)
            {
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
                    if (diff > limit)
                        break;
                }
            }
            return diff;
        }

        public static int CountYCutHighMannose(IGlycan sub, IGlycan glycan, int limit)
        {
            int[] table = glycan.Table();
            int[] subTable = sub.Table();

            // partial core
            int diff = 0;
            if (subTable[0] + subTable[1] + subTable[2] < 3)
            {
                if (subTable[0] == 0)
                    return -1;
                // fucose core
                if (subTable[2] != table[2])
                    diff++;
                diff++;
                return diff;
            }

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
                diff++;

            // branch
            if (subTable[3] > 0)
            {
                for (int i = 0; i < 2; i++)
                {
                    if (table[i + 5] != subTable[i + 5])
                        diff++;
                    if (diff > limit)
                        break;
                }
            }
            if (subTable[4] > 0 && table[7] != subTable[7])
            {
                diff++;
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
            foreach(Monosaccharide sugar in subCompose.Keys)
            {
                compose[sugar] -= subCompose[sugar];
                if (compose[sugar] == 0)
                {
                    compose.Remove(sugar);
                }
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


        public static bool ContainsCut(IGlycan glycan, IGlycan sub, IGlycan subSub)
        {
            int diff1 = CountYCut(sub, glycan, 1);
            int diff2 = CountYCut(subSub, glycan, 1);

            return diff1 == diff2; 
        }

        //// stackoverflow code:
        //// https://stackoverflow.com/questions/756055/listing-all-permutations-of-a-string-integer
        //private static IEnumerable<IEnumerable<T>> GetPermutations<T>(
        //   IEnumerable<T> list, int length)
        //{
        //    if (length == 1) return list.Select(t => new T[] { t });
        //    return GetPermutations(list, length - 1)
        //        .SelectMany(t => list.Where(e => !t.Contains(e)),
        //            (t1, t2) => t1.Concat(new T[] { t2 }));
        //}

        //public static List<IGlycan> PermutationComplex(IGlycan glycan)
        //{
        //    Dictionary<string, IGlycan> results = new Dictionary<string, IGlycan>();
        //    // permuate index (1, 2, 3, 4) -- four branches
        //    IEnumerable<IEnumerable<int>> indexPermuted = GetPermutations(Enumerable.Range(0, 4), 4);
        //    foreach (IEnumerable<int> permute in indexPermuted)
        //    {
        //        IGlycan newGlycan = glycan.Clone();
        //        int[] table = glycan.Table();
        //        int[] subTable = newGlycan.Table();
        //        for (int i = 0; i < 4; i++)
        //        {
        //            int index = permute.ToList()[i];
        //            if (index < 2)
        //            {
        //                subTable[4 + i / 2] = table[4];
        //            }
        //            else
        //            {
        //                subTable[4 + i / 2] = table[5];
        //            }
        //            subTable[i + 6] = table[index + 6];
        //            subTable[i + 10] = table[index + 10];
        //            subTable[i + 14] = table[index + 14];
        //            subTable[i + 18] = table[index + 18];
        //            subTable[i + 22] = table[index + 22];
        //        }
        //        newGlycan.SetTable(subTable);
        //        results[newGlycan.ID()] = newGlycan;
        //    }
        //    return results.Values.ToList();
        //}

        //public static List<IGlycan> PermutationHybrid(IGlycan glycan)
        //{
        //    Dictionary<string, IGlycan> results = new Dictionary<string, IGlycan>();
        //    for (int i = 0; i < 2; i++)
        //    {
        //        for (int j = 0; j < 2; j++)
        //        {
        //            IGlycan newGlycan = glycan.Clone();
        //            int[] table = glycan.Table();
        //            int[] subTable = newGlycan.Table();
        //            subTable[4] = table[4 + i];
        //            subTable[5] = table[5 - i];
        //            for (int k = 0; k < 6; k++)
        //            {
        //                subTable[6 + k * 2] = table[6 + k * 2 + j];
        //                subTable[7 + k * 2] = table[7 + k * 2 - j];
        //            }
        //            newGlycan.SetTable(subTable);
        //            results[newGlycan.ID()] = newGlycan;
        //        }
        //    }
        //    return results.Values.ToList();
        //}

        //// permutate the position of 
        //public static List<IGlycan> Permutation(IGlycan glycan)
        //{
        //    List<IGlycan> results = new List<IGlycan>();
        //    switch (glycan.Type())
        //    {
        //        case GlycanType.NGlycanComplex:
        //            return PermutationComplex(glycan);
        //        case GlycanType.NGlycanHybrid:
        //            return PermutationHybrid(glycan);
        //        case GlycanType.NHighMannose:
        //            results.Add(glycan);
        //            break;
        //    }
        //    return results;
        //}
    }
}
