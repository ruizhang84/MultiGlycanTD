using MultiGlycanTDLibrary.algorithm;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Shapes;

namespace MultiGlycanTD
{
    /// <summary>
    /// Interaction logic for ConfigureWindow.xaml
    /// </summary>
    public partial class ConfigureWindow : Window
    {
        public ConfigureWindow()
        {
            InitializeComponent();
            InitWindow();
        }
        public void InitWindow()
        {
            MS1Tol.Text = SearchingParameters.Access.MS1Tolerance.ToString();
            MSMS2Tol.Text = SearchingParameters.Access.MSMSTolerance.ToString();
            if (SearchingParameters.Access.MS1ToleranceBy
                == ToleranceBy.Dalton)
            {
                MS1TolByDalton.IsChecked = true;
            }
            if (SearchingParameters.Access.MS2ToleranceBy
                == ToleranceBy.Dalton)
            {
                MS2TolByDalton.IsChecked = true;
            }

            foreach (double ion in SearchingParameters.Access.Ions)
            {
                if (ion == MultiGlycanTDLibrary.util.mass.Spectrum.Proton)
                {
                    Proton.IsChecked = true;
                }
                else if (ion == MultiGlycanTDLibrary.util.mass.Spectrum.Sodium)
                {
                    Sodium.IsChecked = true;
                }
            }

            ThreadNums.Text = SearchingParameters.Access.ThreadNums.ToString();

            Cutoff.Text = (SearchingParameters.Access.Cutoff * 100.0).ToString();
            FDR.Text = (SearchingParameters.Access.FDR * 100.0).ToString();
        }

        private void Save_Click(object sender, RoutedEventArgs e)
        {
            if (SaveChanges())
            {
                SearchingParameters.Access.Update();
                Close();
            }
        }

        private void Cancel_Click(object sender, RoutedEventArgs e)
        {
            Close();
        }

        private bool SaveChanges()
        {
            return SaveTolerance() &&
                    SaveOutput() &&
                    SaveSearch();
        }

        private bool SaveSearch()
        {
            if (int.TryParse(ThreadNums.Text, out int nums))
            {
                ConfigureParameters.Access.ThreadNums = nums;
            }
            else
            {
                MessageBox.Show("Thread value is invalid!");
                return false;
            }
            if (Proton.IsChecked == true || Sodium.IsChecked == true)
            {
                ConfigureParameters.Access.Ions.Clear();
                if (Proton.IsChecked == true)
                    ConfigureParameters.Access.Ions
                        .Add(MultiGlycanTDLibrary.util.mass.Spectrum.Proton);
                if (Sodium.IsChecked == true)
                    ConfigureParameters.Access.Ions
                        .Add(MultiGlycanTDLibrary.util.mass.Spectrum.Sodium);
            }
            else
            {
                MessageBox.Show("At least one Ion type!");
                return false;
            }

            return true;
        }


        private bool SaveTolerance()
        {
            if (double.TryParse(MS1Tol.Text, out double tol))
            {
                ConfigureParameters.Access.MS1Tolerance = tol;
            }
            else
            {
                MessageBox.Show("MS tolerance value is invalid!");
                return false;
            }

            if (double.TryParse(MSMS2Tol.Text, out tol))
            {
                ConfigureParameters.Access.MSMSTolerance = tol;
            }
            else
            {
                MessageBox.Show("MSMS tolerance value is invalid!");
                return false;
            }
            return true;
        }

        private bool SaveOutput()
        {
            if (double.TryParse(Cutoff.Text, out double cutoff) && cutoff >= 0 && cutoff <= 100)
            {
                ConfigureParameters.Access.Cutoff = cutoff * 0.01;
            }
            else
            {
                MessageBox.Show("Cutoff level is invalid!");
                return false;
            }

            if (double.TryParse(FDR.Text, out double rate) && rate >= 0 && rate <= 100)
            {
                ConfigureParameters.Access.FDR = rate * 0.01;
            }
            else
            {
                MessageBox.Show("Filter level is invalid!");
                return false;
            }
            return true;
        }

        private void MS1TolByPPM_Checked(object sender, RoutedEventArgs e)
        {
            ConfigureParameters.Access.MS1ToleranceBy = ToleranceBy.PPM;
        }

        private void MS1TolByDalton_Checked(object sender, RoutedEventArgs e)
        {
            ConfigureParameters.Access.MS1ToleranceBy = ToleranceBy.Dalton;
        }

        private void MS2TolByPPM_Checked(object sender, RoutedEventArgs e)
        {
            ConfigureParameters.Access.MS2ToleranceBy = ToleranceBy.PPM;
        }

        private void MS2TolByDalton_Checked(object sender, RoutedEventArgs e)
        {
            ConfigureParameters.Access.MS2ToleranceBy = ToleranceBy.Dalton;
        }
    }
}
