using Microsoft.Win32;
using MultiGlycanTDLibrary.model;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
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

namespace MultiGlycanTD
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
        }

        private void MSMSFileNames_Click(object sender, RoutedEventArgs e)
        {
            OpenFileDialog fileNamesDialog = new OpenFileDialog();
            fileNamesDialog.Filter = "Raw File|*.raw|MGF File|*.mgf";
            fileNamesDialog.Title = "Open a MS2 File";
            fileNamesDialog.Multiselect = true;
            fileNamesDialog.InitialDirectory = Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments);

            if (fileNamesDialog.ShowDialog() == true)
            {
                foreach (string filename in fileNamesDialog.FileNames)
                {
                    if (!SearchingParameters.Access.MSMSFiles.Contains(filename))
                    {
                        lbFiles.Items.Add(filename);
                        SearchingParameters.Access.MSMSFiles.Add(filename);
                    }
                }

            }
        }

        private void DatasetFileNames_Click(object sender, RoutedEventArgs e)
        {
            OpenFileDialog fileNameDialog = new OpenFileDialog();
            fileNameDialog.Filter = "JSON File|*.json";
            fileNameDialog.Title = "Open a JSON File";

            if (fileNameDialog.ShowDialog() == true)
            {
                string jsonStringRead = File.ReadAllText(fileNameDialog.FileName);
                SearchingParameters.Access.Database = JsonSerializer.Deserialize<GlycanJson>(jsonStringRead);
            }
        }
        private void DeselectFiles_Click(object sender, RoutedEventArgs e)
        {
            if (lbFiles.SelectedItem != null)
            {
                string filename = lbFiles.SelectedItem.ToString();
                lbFiles.Items.Remove(lbFiles.SelectedItem);
                if (SearchingParameters.Access.MSMSFiles.Contains(filename))
                    SearchingParameters.Access.MSMSFiles.Remove(filename);
            }
        }

        private void Search_Click(object sender, RoutedEventArgs e)
        {
            if (SearchingParameters.Access.MSMSFiles.Count == 0)
            {
                MessageBox.Show("Please choose MS/MS files");
            }
            else if (SearchingParameters.Access.Database is null)
            {
                MessageBox.Show("Please choose glycan serching (*.JSON) file");
            }
            else
            {
                Window subWindow = new SearchWindow();
                subWindow.Show();
            }
        }

        private void Configure_Click(object sender, RoutedEventArgs e)
        {
            Window subWindow = new ConfigureWindow();
            subWindow.Show();
        }
    }
}
