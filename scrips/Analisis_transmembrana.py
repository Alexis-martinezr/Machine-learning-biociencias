import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import re

class TransmembranePredictor:
    def __init__(self):
        # Escala de hidrofobicidad de Kyte-Doolittle
        self.kyte_doolittle = {
            'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
            'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
            'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
            'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
            'X': 0.0  # Para aminoácidos no estándar
        }
        
        # Escala alternativa más sensible (Wimley-White)
        self.wimley_white = {
            'A': 0.17, 'R': -2.85, 'N': -1.11, 'D': -2.49, 'C': -0.24,
            'Q': -1.10, 'E': -2.68, 'G': -0.01, 'H': -2.33, 'I': 0.31,
            'L': 0.56, 'K': -2.80, 'M': 0.23, 'F': 1.13, 'P': -0.45,
            'S': -0.46, 'T': -0.25, 'W': 1.85, 'Y': 0.94, 'V': -0.07,
            'X': 0.0
        }
        
        # Escala de hidrofobicidad biológica (Hessa et al. 2005) - invertida para hidrofobicidad
        self.hessa_scale = {
            'A': -0.11, 'R': -2.58, 'N': -0.63, 'D': -3.64, 'C': 0.13,
            'Q': -0.77, 'E': -4.16, 'G': -0.74, 'H': -2.06, 'I': 0.60,
            'L': 0.55, 'K': -2.71, 'M': 0.10, 'F': 0.32, 'P': -2.23,
            'S': -0.13, 'T': -0.14, 'W': 0.24, 'Y': -0.68, 'V': 0.31,
            'X': 0.0
        }
        
        # Configuraciones preestablecidas
        self.presets = {
            'conservative': {'threshold': 1.5, 'min_length': 18, 'window': 19, 'scale': 'kyte_doolittle'},
            'standard': {'threshold': 1.0, 'min_length': 15, 'window': 19, 'scale': 'kyte_doolittle'},
            'sensitive': {'threshold': 0.6, 'min_length': 13, 'window': 17, 'scale': 'kyte_doolittle'},
            'ultra_sensitive': {'threshold': 0.3, 'min_length': 11, 'window': 15, 'scale': 'combined'}
        }
        
    def parse_pdb(self, pdb_content):
        """
        Extrae la secuencia de aminoácidos de un archivo PDB
        """
        sequences = defaultdict(list)
        
        for line in pdb_content.split('\n'):
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                chain_id = line[21]
                residue_num = int(line[22:26].strip())
                residue_name = line[17:20].strip()
                
                # Convertir código de 3 letras a 1 letra
                aa_code = self.three_to_one(residue_name)
                sequences[chain_id].append((residue_num, aa_code))
        
        # Ordenar por número de residuo y crear secuencias continuas
        final_sequences = {}
        for chain_id, residues in sequences.items():
            residues.sort(key=lambda x: x[0])
            sequence = ''.join([aa for _, aa in residues])
            final_sequences[chain_id] = sequence
            
        return final_sequences
    
    def three_to_one(self, three_letter_code):
        """
        Convierte código de aminoácido de 3 letras a 1 letra
        """
        conversion = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
        }
        return conversion.get(three_letter_code, 'X')
    
    def calculate_hydropathy(self, sequence, scale='kyte_doolittle', window_size=19):
        """
        Calcula el perfil de hidrofobicidad usando ventana deslizante
        """
        if len(sequence) < window_size:
            return []
        
        # Seleccionar escala
        if scale == 'wimley_white':
            hydro_scale = self.wimley_white
        elif scale == 'hessa':
            hydro_scale = self.hessa_scale
        elif scale == 'combined':
            # Combinar múltiples escalas
            return self.calculate_combined_hydropathy(sequence, window_size)
        else:
            hydro_scale = self.kyte_doolittle
        
        hydropathy_scores = []
        half_window = window_size // 2
        
        for i in range(len(sequence)):
            start = max(0, i - half_window)
            end = min(len(sequence), i + half_window + 1)
            window_seq = sequence[start:end]
            
            # Calcular promedio de hidrofobicidad en la ventana
            scores = [hydro_scale.get(aa, 0.0) for aa in window_seq]
            avg_score = sum(scores) / len(scores) if scores else 0.0
            hydropathy_scores.append(avg_score)
            
        return hydropathy_scores
    
    def calculate_combined_hydropathy(self, sequence, window_size):
        """
        Calcula hidrofobicidad combinando múltiples escalas
        """
        hydro_kd = self.calculate_hydropathy(sequence, 'kyte_doolittle', window_size)
        hydro_ww = self.calculate_hydropathy(sequence, 'wimley_white', window_size)
        hydro_hessa = self.calculate_hydropathy(sequence, 'hessa', window_size)
        
        # Promedio ponderado: KD tiene más peso, luego WW, luego Hessa
        combined = []
        for i in range(len(hydro_kd)):
            combined_score = (2 * hydro_kd[i] + hydro_ww[i] + 0.5 * hydro_hessa[i]) / 3.5
            combined.append(combined_score)
        
        return combined
    
    def predict_transmembrane_regions(self, sequence, mode='standard'):
        """
        Predice regiones transmembranales con diferentes modos de sensibilidad
        """
        if mode not in self.presets:
            mode = 'standard'
            
        preset = self.presets[mode]
        threshold = preset['threshold']
        min_length = preset['min_length']
        window_size = preset['window']
        scale = preset['scale']
        
        # Calcular hidrofobicidad
        hydropathy_scores = self.calculate_hydropathy(sequence, scale, window_size)
        
        if not hydropathy_scores:
            return [], hydropathy_scores
        
        # Suavizado adaptativo
        smoothed_scores = self.smooth_scores(hydropathy_scores, self.get_smoothing_window(mode))
        
        # Encontrar regiones TM
        tm_regions = self.find_tm_regions_advanced(smoothed_scores, threshold, min_length)
        
        # Post-procesamiento según el modo
        if mode in ['sensitive', 'ultra_sensitive']:
            tm_regions = self.merge_close_regions(tm_regions, max_gap=8)
            tm_regions = self.refine_weak_regions(sequence, tm_regions, smoothed_scores, threshold * 0.8)
        else:
            tm_regions = self.merge_close_regions(tm_regions, max_gap=10)
        
        return tm_regions, hydropathy_scores
    
    def get_smoothing_window(self, mode):
        """
        Retorna tamaño de ventana de suavizado según el modo
        """
        smoothing_windows = {
            'conservative': 5,
            'standard': 3,
            'sensitive': 3,
            'ultra_sensitive': 2
        }
        return smoothing_windows.get(mode, 3)
    
    def smooth_scores(self, scores, window=3):
        """
        Suaviza los scores usando media móvil
        """
        if len(scores) < window:
            return scores
            
        smoothed = []
        half_win = window // 2
        
        for i in range(len(scores)):
            start = max(0, i - half_win)
            end = min(len(scores), i + half_win + 1)
            avg = sum(scores[start:end]) / (end - start)
            smoothed.append(avg)
        
        return smoothed
    
    def find_tm_regions_advanced(self, scores, threshold, min_length):
        """
        Método avanzado para encontrar regiones TM
        """
        tm_regions = []
        
        # Identificar todos los puntos sobre el umbral
        above_threshold = [i for i, score in enumerate(scores) if score >= threshold]
        
        if not above_threshold:
            return tm_regions
        
        # Agrupar puntos consecutivos o casi consecutivos
        current_region = [above_threshold[0]]
        max_gap = 5  # Máximo gap permitido
        
        for i in range(1, len(above_threshold)):
            if above_threshold[i] - above_threshold[i-1] <= max_gap:
                current_region.append(above_threshold[i])
            else:
                # Procesar región actual
                if len(current_region) >= min_length:
                    start = current_region[0] + 1  # +1 para numeración de residuos
                    end = current_region[-1] + 1
                    tm_regions.append((start, end))
                current_region = [above_threshold[i]]
        
        # No olvidar la última región
        if len(current_region) >= min_length:
            start = current_region[0] + 1
            end = current_region[-1] + 1
            tm_regions.append((start, end))
        
        return tm_regions
    
    def merge_close_regions(self, regions, max_gap=10):
        """
        Fusiona regiones TM que están muy cerca entre sí
        """
        if len(regions) <= 1:
            return regions
        
        merged = [regions[0]]
        
        for current_start, current_end in regions[1:]:
            last_start, last_end = merged[-1]
            
            # Si la región actual está cerca de la anterior, fusionar
            if current_start - last_end <= max_gap:
                merged[-1] = (last_start, current_end)
            else:
                merged.append((current_start, current_end))
        
        return merged
    
    def refine_weak_regions(self, sequence, tm_regions, scores, weak_threshold):
        """
        Refina regiones débiles verificando composición de aminoácidos
        """
        refined_regions = []
        
        for start, end in tm_regions:
            # Extraer secuencia de la región
            region_seq = sequence[start-1:end]
            
            # Calcular porcentaje de aminoácidos hidrofóbicos
            hydrophobic_aa = {'A', 'V', 'L', 'I', 'M', 'F', 'W', 'Y', 'C'}
            hydrophobic_percent = sum(1 for aa in region_seq if aa in hydrophobic_aa) / len(region_seq)
            
            # Calcular score promedio de la región
            region_scores = scores[start-1:end]
            avg_score = sum(region_scores) / len(region_scores) if region_scores else 0
            
            # Criterios para mantener la región
            keep_region = (
                avg_score >= weak_threshold or      # Score promedio sobre umbral débil
                hydrophobic_percent >= 0.5 or      # Al menos 50% hidrofóbicos
                len(region_seq) >= 20               # Longitud larga compensa score bajo
            )
            
            if keep_region:
                refined_regions.append((start, end))
        
        return refined_regions
    
    def analyze_protein(self, pdb_content, modes=['standard']):
        """
        Análisis completo de una proteína con múltiples modos
        """
        sequences = self.parse_pdb(pdb_content)
        results = {}
        
        for chain_id, sequence in sequences.items():
            if len(sequence) > 0:
                chain_results = {
                    'sequence': sequence,
                    'length': len(sequence),
                    'hydrophobic_content': self.calculate_hydrophobic_content(sequence)
                }
                
                # Análisis con cada modo solicitado
                for mode in modes:
                    tm_regions, hydropathy_scores = self.predict_transmembrane_regions(sequence, mode)
                    chain_results[f'tm_regions_{mode}'] = tm_regions
                    chain_results[f'num_tm_helices_{mode}'] = len(tm_regions)
                    
                    if mode == modes[0]:  # Guardar scores del primer modo para graficar
                        chain_results['hydropathy_scores'] = hydropathy_scores
                
                results[chain_id] = chain_results
        
        return results
    
    def calculate_hydrophobic_content(self, sequence):
        """
        Calcula el porcentaje de aminoácidos hidrofóbicos
        """
        hydrophobic_aa = {'A', 'V', 'L', 'I', 'M', 'F', 'W', 'Y', 'C'}
        hydrophobic_count = sum(1 for aa in sequence if aa in hydrophobic_aa)
        return (hydrophobic_count / len(sequence)) * 100 if sequence else 0
    
    def compare_all_modes(self, pdb_content):
        """
        Compara todos los modos de predicción
        """
        all_modes = ['conservative', 'standard', 'sensitive', 'ultra_sensitive']
        return self.analyze_protein(pdb_content, all_modes)
    
    def plot_hydropathy_multi_mode(self, sequence, hydropathy_scores, all_tm_regions, chain_id, modes):
        """
        Visualiza el perfil de hidrofobicidad con múltiples modos
        """
        plt.figure(figsize=(16, 10))
        
        # Gráfico de hidrofobicidad
        positions = list(range(1, len(hydropathy_scores) + 1))
        plt.plot(positions, hydropathy_scores, 'b-', linewidth=1.5, alpha=0.8, label='Perfil de hidrofobicidad')
        
        # Líneas de umbral para cada modo
        colors = {'conservative': 'red', 'standard': 'orange', 'sensitive': 'green', 'ultra_sensitive': 'purple'}
        for mode in modes:
            if mode in self.presets:
                threshold = self.presets[mode]['threshold']
                plt.axhline(y=threshold, color=colors[mode], linestyle='--', alpha=0.7, 
                           label=f'Umbral {mode} ({threshold})')
        
        plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
        
        # Destacar regiones transmembranales por modo
        y_offset = 0
        for mode in modes:
            if f'tm_regions_{mode}' in all_tm_regions:
                regions = all_tm_regions[f'tm_regions_{mode}']
                for i, (start, end) in enumerate(regions):
                    alpha = 0.3 - y_offset * 0.05
                    plt.axvspan(start, end, alpha=max(alpha, 0.1), color=colors[mode], 
                               label=f'TMD {mode}' if i == 0 else "")
                y_offset += 1
        
        plt.xlabel('Posición del Residuo')
        plt.ylabel('Hidrofobicidad (Kyte-Doolittle)')
        plt.title(f'Predicción Comparativa de Regiones Transmembranales - Cadena {chain_id}')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True, alpha=0.3)
        
        # Añadir información resumida
        info_text = "TMDs detectados:\n"
        for mode in modes:
            if f'num_tm_helices_{mode}' in all_tm_regions:
                count = all_tm_regions[f'num_tm_helices_{mode}']
                info_text += f"{mode}: {count}\n"
        
        plt.text(0.02, 0.98, info_text, transform=plt.gca().transAxes, 
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        plt.show()
    
    def print_results(self, results):
        """
        Imprime resultados del análisis con todos los modos
        """
        print("=" * 80)
        print("ANÁLISIS COMPARATIVO DE REGIONES TRANSMEMBRANALES")
        print("=" * 80)
        
        for chain_id, data in results.items():
            print(f"\n Cadena {chain_id}:")
            print(f"   Longitud: {data['length']} residuos")
            print(f"   Hidrofobicidad: {data['hydrophobic_content']:.1f}%")
            
            # Mostrar resultados por modo
            modes = ['conservative', 'standard', 'sensitive', 'ultra_sensitive']
            mode_names = {'conservative': 'Conservador', 'standard': 'Estándar', 
                         'sensitive': 'Sensible', 'ultra_sensitive': 'Ultra-sensible'}
            
            print(f"\n   Resumen de TMDs detectados:")
            for mode in modes:
                key = f'num_tm_helices_{mode}'
                if key in data:
                    count = data[key]
                    print(f"      {mode_names[mode]:12}: {count} TMDs")
            
            # Mostrar detalles de cada modo con TMDs
            for mode in modes:
                regions_key = f'tm_regions_{mode}'
                if regions_key in data and data[regions_key]:
                    print(f"\n   Regiones {mode_names[mode].lower()}:")
                    for i, (start, end) in enumerate(data[regions_key], 1):
                        length = end - start + 1
                        tm_seq = data['sequence'][start-1:end]
                        
                        # Calcular métricas de la región
                        hydrophobic_aa = {'A', 'V', 'L', 'I', 'M', 'F', 'W', 'Y', 'C'}
                        hydrophobic_percent = sum(1 for aa in tm_seq if aa in hydrophobic_aa) / len(tm_seq) * 100
                        avg_hydro = sum(self.kyte_doolittle.get(aa, 0) for aa in tm_seq) / len(tm_seq)
                        
                        print(f"      TM{i}: {start:3d}-{end:3d} ({length:2d} aa) "
                              f"Hidro:{hydrophobic_percent:4.1f}% Score:{avg_hydro:5.2f}")
                        print(f"           {tm_seq}")
        
        print("\n" + "=" * 80)
        print("   Recomendación: Para proteínas con 8 TMDs conocidos,")
        print("   usa el modo 'ultra_sensitive' o ajusta manualmente los parámetros.")
        print("=" * 80)

# Funciones de análisis
def analyze_pdb_comprehensive(pdb_file_path):
    """
    Análisis completo con todos los modos para encontrar el óptimo
    """
    predictor = TransmembranePredictor()
    
    try:
        with open(pdb_file_path, 'r') as file:
            pdb_content = file.read()
        
        # Análisis con todos los modos
        results = predictor.compare_all_modes(pdb_content)
        predictor.print_results(results)
        
        # Generar gráfico comparativo
        for chain_id, data in results.items():
            modes = ['conservative', 'standard', 'sensitive', 'ultra_sensitive']
            predictor.plot_hydropathy_multi_mode(
                data['sequence'], 
                data['hydropathy_scores'], 
                data, 
                chain_id, 
                modes
            )
        
        return results
        
    except FileNotFoundError:
        print(f" Error: No se pudo encontrar el archivo {pdb_file_path}")
        return None
    except Exception as e:
        print(f" Error al procesar el archivo: {str(e)}")
        return None

def analyze_pdb_file(pdb_file_path, mode='standard'):
    """
    Análisis con modo específico
    """
    predictor = TransmembranePredictor()
    
    try:
        with open(pdb_file_path, 'r') as file:
            pdb_content = file.read()
        
        results = predictor.analyze_protein(pdb_content, [mode])
        predictor.print_results(results)
        
        # Graficar
        for chain_id, data in results.items():
            if data['hydropathy_scores']:
                predictor.plot_hydropathy_multi_mode(
                    data['sequence'], 
                    data['hydropathy_scores'], 
                    data, 
                    chain_id,
                    [mode]
                )
        
        return results
        
    except FileNotFoundError:
        print(f" Error: No se pudo encontrar el archivo {pdb_file_path}")
        return None
    except Exception as e:
        print(f" Error al procesar el archivo: {str(e)}")
        return None

# Ejemplo de uso con secuencia directa
def example_with_sequence():
    """
    Ejemplo usando una secuencia conocida
    """
    # Secuencia ejemplo de una proteína transmembranal
    example_sequence = "MKTFFLVLLLLVLPVLSSVHHHHHHENLYFQGAMGSEFINDRIYQNADAYVDSLKGPFLLAIILVAVGVFFVGALFLEYRRKLRKRWQIKGTLDQIRAGLEFVDKEQPFSAQILASLGGPKILLILSALSMLVFGVVFYEHDCNMVDIFMVLISAVVFLSQSFLKLYRFGDAVCSLLFVFFMGFALPLLILGDWLALVFQVPLAVVFVAFFFGWLYRWNSHWKLPSLRRYLPFGRGLLVLVGAALAMLALYACQPLAGSLIIYFMLLGWLFYQGKKLILGLWVIIFSWVVGILPVVVFVCQSGNNLMTQSLVLPHLQKGFWRTSALLFFVWMEGQVLRFLNQASLWTLRGILMLVFFGGGSQRYRFFYILQGTGVSRRQTFGVPGPALFVWWYIGRTVIAVPFFFQGFTFLQSFTGFFMSLFQLLGAFVLLAFYLQKQVNIYQKLIRQVCKGKEELARQHVQAMQ"
    
    predictor = TransmembranePredictor()
    
    print("Ejemplo con secuencia larga (posible proteína multi-TMD):")
    print(f"Longitud: {len(example_sequence)} residuos")
    
    # Probar todos los modos
    modes = ['conservative', 'standard', 'sensitive', 'ultra_sensitive']
    for mode in modes:
        tm_regions, hydropathy_scores = predictor.predict_transmembrane_regions(example_sequence, mode)
        print(f"\nModo {mode}: {len(tm_regions)} TMDs detectados")
        for i, (start, end) in enumerate(tm_regions, 1):
            print(f"  TM{i}: {start}-{end} ({end-start+1} aa)")

if __name__ == "__main__":
    print("🧬 Predictor Mejorado de Regiones Transmembranales")
    print("=" * 60)
    print("Funciones disponibles:")
    print("1. analyze_pdb_file('archivo.pdb', 'mode') - Análisis específico")
    print("   Modos: 'conservative', 'standard', 'sensitive', 'ultra_sensitive'")
    print("2. analyze_pdb_comprehensive('archivo.pdb') - Comparar todos los modos")
    print("3. example_with_sequence() - Ejemplo con secuencia")
    print("\n💡 Para detectar 8 TMDs, prueba 'ultra_sensitive' o 'sensitive'")
    print("=" * 60)
    
    # Ejecutar ejemplo
    example_with_sequence()

import pandas as pd

def results_to_dataframe(results):
    """Convierte results a DataFrame para fácil exploración"""
    rows = []
    
    for chain_id, data in results.items():
        # Información básica
        base_row = {
            'chain_id': chain_id,
            'sequence': data['sequence'],
            'length': data['length'],
            'hydrophobic_content': data['hydrophobic_content']
        }
        
        # Agregar TMDs por modo
        modes = ['conservative', 'standard', 'sensitive', 'ultra_sensitive']
        for mode in modes:
            tm_key = f'tm_regions_{mode}'
            num_key = f'num_tm_helices_{mode}'
            if tm_key in data:
                base_row[f'tmds_{mode}'] = data[num_key]
                base_row[f'regions_{mode}'] = str(data[tm_key])
        
        rows.append(base_row)
    
    return pd.DataFrame(rows)

from Bio.PDB import PDBParser

def structur_protein_df(pdb_file):

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    data = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id()[0] == " ":  # Aminoácidos estándar
                    if residue.has_id("CA"):
                        ca = residue["CA"]
                        data.append({
                            "model": model.id,
                            "chain": chain.id,
                            "residue_name": residue.get_resname(),
                            "residue_id": residue.get_id()[1],
                            "x": ca.coord[0],
                            "y": ca.coord[1],
                            "z": ca.coord[2]
                        })

    return pd.DataFrame(data)

import ast

def aa(transmembrana, cordenadas_proteina):
    # Supongamos que la columna tiene una sola fila con esa string
    tmds = [list(t) for t in ast.literal_eval(transmembrana['regions_ultra_sensitive'].iloc[0])]

    
    aa =[]
    for a in tmds:
        rango = list(range(a[0], a[1]+1))
        aa.extend(rango)


    region = []
    for a in cordenadas_proteina['residue_id']:
        if a in aa:
            region.append('Transmembrana')
        else:
            region.append('No_transmembrana')
    return region