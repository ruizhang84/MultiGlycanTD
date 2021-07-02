using Microsoft.Win32;
using MultiGlycanClassLibrary.util.mass;
using MultiGlycanTDLibrary.engine.glycan;
using MultiGlycanTDLibrary.model;
using MultiGlycanTDLibrary.model.glycan;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace GlycanCalculator
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        protected string fileName;
        protected string glycanFileName = "";
        protected int hexNAc;
        protected int hex;
        protected int fuc;
        protected int neuAc;
        protected int neuGc;
        protected bool complexInclude;
        protected bool hybridInclude;
        protected bool highMannoseInclude;
        protected int order;
        protected int precision;

        protected Derivatization derivatization;
        protected bool permethylated;
        protected bool reduced;
        protected List<FragmentTypes> types = new List<FragmentTypes>();

        public MainWindow()
        {
            InitializeComponent();
        }

        private void JSONFileNames_Click(object sender, RoutedEventArgs e)
        {
            SaveFileDialog savefile = new SaveFileDialog();
            savefile.FileName = "database.json";
            savefile.Filter = "Json files (*.json) | *.json";

            if (savefile.ShowDialog() == true)
            {
                fileName = savefile.FileName;
            }
        }

        private void TxtFileNames_Click(object sender, RoutedEventArgs e)
        {
            OpenFileDialog openfile = new OpenFileDialog();
            openfile.Filter = "Text files (*.txt) | *.txt";

            if (openfile.ShowDialog() == true)
            {
                glycanFileName = openfile.FileName;
            }
        }

        private bool SearchInit()
        {
            int bound = 12;
            if (int.TryParse(HexNAc.Text, out bound) && bound >= 0)
            {
                hexNAc = bound;
            }
            else
            {
                MessageBox.Show("HexNAc value is invalid!");
                return false;
            }
            if (int.TryParse(Hex.Text, out bound) && bound >= 0)
            {
                hex = bound;
            }
            else
            {
                MessageBox.Show("HexNAc value is invalid!");
                return false;
            }
            if (int.TryParse(Fuc.Text, out bound) && bound >= 0)
            {
                fuc = bound;
            }
            else
            {
                MessageBox.Show("Fuc value is invalid!");
                return false;
            }
            if (int.TryParse(NeuAc.Text, out bound) && bound >= 0)
            {
                neuAc = bound;
            }
            else
            {
                MessageBox.Show("NeuAc value is invalid!");
                return false;
            }
            if (int.TryParse(NeuGc.Text, out bound) && bound >= 0)
            {
                neuGc = bound;
            }
            else
            {
                MessageBox.Show("NeuGc value is invalid!");
                return false;
            }
            if (int.TryParse(Order.Text, out bound) && bound >= 0)
            {
                order = bound;
            }
            else
            {
                MessageBox.Show("Num Distribution value is invalid!");
                return false;
            }
            if (int.TryParse(Precision.Text, out bound) && bound >= 0)
            {
                precision = bound;
            }
            else
            {
                MessageBox.Show("Decimal value is invalid!");
                return false;
            }
            if (ComplexNGlycan.IsChecked == false &&
                HybridNGlycan.IsChecked == false && HighMannose.IsChecked == false)
            {
                MessageBox.Show("Choose at least one NGlycan type!");
                return false;
            }
            else
            {
                complexInclude = ComplexNGlycan.IsChecked == true;
                hybridInclude = HybridNGlycan.IsChecked == true;
                highMannoseInclude = HighMannose.IsChecked == true;
            }

            types.Clear();
            if (Bions.IsChecked == true)
                types.Add(FragmentTypes.B);
            if (Cions.IsChecked == true)
                types.Add(FragmentTypes.C);
            if (Yions.IsChecked == true)
                types.Add(FragmentTypes.Y);
            if (Zions.IsChecked == true)
                types.Add(FragmentTypes.Z);
            if (BYions.IsChecked == true)
                types.Add(FragmentTypes.BY);
            if (BZions.IsChecked == true)
                types.Add(FragmentTypes.BZ);
            if (CYions.IsChecked == true)
                types.Add(FragmentTypes.CY);
            if (YYions.IsChecked == true)
                types.Add(FragmentTypes.YY);
            if (YZions.IsChecked == true)
                types.Add(FragmentTypes.YZ);
            if (ZZions.IsChecked == true)
                types.Add(FragmentTypes.ZZ); 
            if (types.Count == 0)
            {
                MessageBox.Show("Select at least one ions!");
                return false;
            }
            
            if (fileName is null || fileName.Length == 0)
            {
                MessageBox.Show("Name the saving file!");
                return false;
            }

            permethylated = Permethylated.IsChecked == true;
            reduced = PermethylatedReduced.IsChecked == true;
            if (unDerived.IsChecked == true)
            {
                derivatization = Derivatization.Underivatized;
            }
            else if (o2AA.IsChecked == true)
            {
                derivatization = Derivatization.k2AA;
            }
            else if (o2AB.IsChecked == true)
            {
                derivatization = Derivatization.k2AB;
                
            }
            return true;
        }

        private Task Process()
        {
            // init mass calculator
            if (permethylated)
            {
                GlycanIonsBuilder.Build.Permethylated = true;
                Glycan.To.SetPermethylation(true, reduced);
            }
            else
            {
                GlycanIonsBuilder.Build.Permethylated = false;
                Glycan.To.SetPermethylation(false, reduced);
                switch (derivatization)
                {
                    case Derivatization.k2AA:
                        GlycanIonsBuilder.Build.Derivatization = GlycanIonsBuilder.k2AA;
                        Glycan.To.Derivatization = Glycan.k2AA;
                        break;
                    case Derivatization.k2AB:
                        GlycanIonsBuilder.Build.Derivatization = GlycanIonsBuilder.k2AB;
                        Glycan.To.Derivatization = Glycan.k2AB;
                        break;
                    default:
                        GlycanIonsBuilder.Build.Derivatization = GlycanIonsBuilder.kWater;
                        Glycan.To.Derivatization = Glycan.kWater;
                        break;
                }
            }
            GlycanIonsBuilder.Build.Types = types;

            // build
            IGlycanBuilder glycanBuilder;
            if (glycanFileName.Length > 0)
            {
                List<SortedDictionary<Monosaccharide, int>> glycanList
                    =  CalculatorHelper.ReadFilter(glycanFileName);

                glycanBuilder =
                new GlycanBuilderFiltered(glycanList, hexNAc, hex, fuc, neuAc, neuGc,
                complexInclude, hybridInclude, highMannoseInclude,
                order, permethylated, reduced, derivatization);
            }
            else
            {
                glycanBuilder =
                new GlycanBuilder(hexNAc, hex, fuc, neuAc, neuGc,
                complexInclude, hybridInclude, highMannoseInclude,
                order, permethylated, reduced, derivatization);
            }

            glycanBuilder.Build();

            // distribution maps
            var distr_map = glycanBuilder.GlycanDistribMaps();
            var mass_map = glycanBuilder.GlycanMassMaps();
            CompdJson compdJson = new CompdJson()
            {
                DistrMap = distr_map,
                MassMap = mass_map
            };

            // fragmentation maps
            object obj = new object();
            Dictionary<double, Dictionary<FragmentTypes, List<string>>> fragments
                = new Dictionary<double, Dictionary<FragmentTypes, List<string>>>();
            List<Tuple<string, FragmentTypes, List<double>>> fragmentsContainer
                = new List<Tuple<string, FragmentTypes, List<double>>>();
            var map = glycanBuilder.GlycanMaps();


            Parallel.ForEach(map, pair =>
            {
                var id = pair.Key;
                var glycan = pair.Value;
                if (glycan.IsValid())
                {
                    List<Tuple<string, FragmentTypes, List<double>>> fragmentMass
                        = new List<Tuple<string, FragmentTypes, List<double>>>();
                    foreach (FragmentTypes type in GlycanIonsBuilder.Build.Types)
                    {
                        List<double> massList = GlycanIonsBuilder.Build.Fragments(glycan, type)
                                        .Select(m => Math.Round(m, 4)).ToList();
                        fragmentMass.Add(Tuple.Create(id, type, massList));
                    }

                    lock (obj)
                    {
                        fragmentsContainer.AddRange(fragmentMass);
                    }
                }
            });

            foreach (Tuple<string, FragmentTypes, List<double>> item in fragmentsContainer)
            {
                string id = item.Item1;
                FragmentTypes type = item.Item2;

                foreach (double mass in item.Item3)
                {
                    if (!fragments.ContainsKey(mass))
                        fragments[mass] = new Dictionary<FragmentTypes, List<string>>();
                    if (!fragments[mass].ContainsKey(type))
                        fragments[mass][type] = new List<string>();
                    fragments[mass][type].Add(id);
                }
            }

            fragmentsContainer.Clear();
            var composition_map = glycanBuilder.GlycanCompositionMaps();
            Dictionary<string, List<string>> id_map = new Dictionary<string, List<string>>();
            foreach (var pair in composition_map)
            {
                id_map[pair.Key] = pair.Value.Select(p => p.ID()).ToList();
            }

            // serialization
            DerivationType derivation = DerivationType.Native;
            if (permethylated)
            {
                derivation = DerivationType.Permethylated;
            }

            GlycanJson glycanJson = new GlycanJson()
            {
                Derivation = derivation,
                Compound = compdJson,
                IDMap = id_map,
                FragmentMap = fragments
            };
            string jsonString = JsonSerializer.Serialize(glycanJson);
            File.WriteAllText(fileName, jsonString);

            return Task.CompletedTask;
        }

        private void Search_Click(object sender, RoutedEventArgs e)
        {
            if (!SearchInit())
                return;

            Search.IsEnabled = false;
            Task.Run(Process);
        }

        private void Permethylated_Checked(object sender, RoutedEventArgs e)
        {
            permethylated = true;
            if (PermethylatedReduced.IsEnabled == false)
            {
                PermethylatedReduced.IsEnabled = true;
            }
            if (NativeDerivatization.IsEnabled == true)
            {
                NativeDerivatization.IsEnabled = false;
            }
        }

        private void Native_Checked(object sender, RoutedEventArgs e)
        {
            permethylated = false;
            if (PermethylatedReduced.IsEnabled == true)
            {
                PermethylatedReduced.IsEnabled = false;
            }
            if (NativeDerivatization.IsEnabled == false)
            {
                NativeDerivatization.IsEnabled = true;
            }
        }

        private void unDerived_Checked(object sender, RoutedEventArgs e)
        {
            derivatization = Derivatization.Underivatized;
        }

        private void o2AA_Checked(object sender, RoutedEventArgs e)
        {
            derivatization = Derivatization.k2AA;
        }

        private void o2AB_Checked(object sender, RoutedEventArgs e)
        {
            derivatization = Derivatization.k2AB;
        }

    }
}
