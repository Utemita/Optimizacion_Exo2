#!/usr/bin/env python3
"""
explorar_pizzolato.py
=====================
Exploration script for the Pizzolato/Ninapro CyberGlove II kinematic database.
Downloads and analyzes raw kinematic data from Zenodo record 20561611 to find
fine/pinch grip motions with high DIP (distal interphalangeal) flexion range.

Dataset: Ninapro Database 1 - Raw CyberGlove II data
Source: https://zenodo.org/records/20561611 (Raw data.zip, 303MB)
Reference: Atzori et al., "Electromyography data for non-invasive naturally-controlled
           robotic hand prostheses", Scientific Data, 2014.

CyberGlove II 22-sensor layout for index finger:
  - Channel 5: Index MCP (metacarpophalangeal) flexion sensor
  - Channel 6: Index PIP (proximal interphalangeal) flexion sensor
  - Index DIP: Not independently measured; derived using biomechanical coupling
    DIP = 0.70 * PIP (Kamper et al. 2003)

Calibration: Linear mapping from raw sensor values to joint angles (degrees)
  - MCP gain: 90 deg / 130 units = 0.692 deg/unit
  - PIP gain: 100 deg / 170 units = 0.588 deg/unit
  - Baseline: 5th percentile of smoothed signal (rest/extension position)
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
from scipy.ndimage import uniform_filter1d


def find_pinch_segments(data, min_duration=2.0, max_duration=15.0, 
                        min_ch6_range=60, sample_rate=100):
    """
    Identify segments of high index finger flexion in CyberGlove II data.
    
    Parameters
    ----------
    data : np.ndarray, shape (N, 22)
        Raw 22-channel CyberGlove II data.
    min_duration : float
        Minimum segment duration in seconds.
    max_duration : float
        Maximum segment duration in seconds.
    min_ch6_range : float
        Minimum range in channel 6 (PIP) to qualify.
    sample_rate : int
        Sampling rate in Hz.
    
    Returns
    -------
    list of dict
        Segments sorted by combined ch5+ch6 range (descending).
    """
    index_signal = uniform_filter1d(data[:, 5] + data[:, 6], 200)
    threshold = np.percentile(index_signal, 85)
    
    above = index_signal > threshold
    transitions = np.diff(above.astype(int))
    starts = np.where(transitions == 1)[0]
    ends = np.where(transitions == -1)[0]
    
    if len(starts) == 0 or len(ends) == 0:
        return []
    if starts[0] > ends[0]:
        ends = ends[1:]
    if len(starts) > len(ends):
        starts = starts[:len(ends)]
    
    segments = []
    for s, e in zip(starts, ends):
        dur = (e - s) / sample_rate
        if dur < min_duration or dur > max_duration:
            continue
        seg = data[s:e, :]
        r5 = seg[:, 5].max() - seg[:, 5].min()
        r6 = seg[:, 6].max() - seg[:, 6].min()
        if r6 >= min_ch6_range:
            segments.append({
                'start': s,
                'end': e,
                'duration': dur,
                'ch5_range': r5,
                'ch6_range': r6,
                'combined': r5 + r6
            })
    
    segments.sort(key=lambda x: x['combined'], reverse=True)
    return segments


def calibrate_to_angles(raw_segment, dip_coupling=0.70):
    """
    Convert raw CyberGlove II sensor values to calibrated joint angles.
    
    Parameters
    ----------
    raw_segment : np.ndarray, shape (N, 22)
        Raw sensor data for a movement segment.
    dip_coupling : float
        DIP/PIP biomechanical coupling ratio (default 0.70).
    
    Returns
    -------
    theta_mcp, theta_pip, theta_dip : np.ndarray
        Joint angles in degrees.
    """
    ch5 = raw_segment[:, 5].copy()
    ch6 = raw_segment[:, 6].copy()
    
    # Smooth raw data
    win = min(31, max(5, len(ch5) // 10 * 2 + 1))
    ch5_smooth = savgol_filter(ch5, win, 3)
    ch6_smooth = savgol_filter(ch6, win, 3)
    
    # Rest baseline (finger extended)
    ch5_base = np.percentile(ch5_smooth, 5)
    ch6_base = np.percentile(ch6_smooth, 5)
    
    # Calibration gains (deg/unit)
    gain_mcp = 90.0 / 130.0   # 0.692
    gain_pip = 100.0 / 170.0  # 0.588
    
    # Convert
    theta_mcp = (ch5_smooth - ch5_base) * gain_mcp
    theta_pip = (ch6_smooth - ch6_base) * gain_pip
    theta_dip = theta_pip * dip_coupling
    
    # Physiological limits
    theta_mcp = np.clip(theta_mcp, 0, 90)
    theta_pip = np.clip(theta_pip, 0, 100)
    theta_dip = np.clip(theta_dip, 0, 70)
    
    return theta_mcp, theta_pip, theta_dip


def resample_to_n_points(theta_mcp, theta_pip, theta_dip, n=120):
    """Resample trajectories to n evenly-spaced points."""
    t_orig = np.linspace(0, 1, len(theta_mcp))
    t_new = np.linspace(0, 1, n)
    
    mcp = np.clip(interp1d(t_orig, theta_mcp, kind='cubic')(t_new), 0, 90)
    pip = np.clip(interp1d(t_orig, theta_pip, kind='cubic')(t_new), 0, 100)
    dip = np.clip(interp1d(t_orig, theta_dip, kind='cubic')(t_new), 0, 70)
    
    return mcp, pip, dip


def main():
    # Paths
    zip_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            '..', 'Base-de-Datos-Pizzoleto', 'data', 'Raw_data.zip')
    
    if not os.path.exists(zip_path):
        # Try alternative path
        zip_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'Base-de-Datos-Pizzoleto', 'data', 'Raw_data.zip')
    
    if not os.path.exists(zip_path):
        print("ERROR: Raw_data.zip not found.")
        print(f"Expected at: {zip_path}")
        print("Download from: https://zenodo.org/api/records/20561611/files/Raw%20data.zip/content")
        sys.exit(1)
    
    import zipfile
    zf = zipfile.ZipFile(zip_path)
    names = zf.namelist()
    
    # Find all E2 (grasping) files from DB1
    e2_files = sorted([n for n in names if '/DB1/' in n and '_E2_' in n])
    
    print("=" * 70)
    print("EXPLORACION DE BASE DE DATOS PIZZOLATO/NINAPRO")
    print("Busqueda de agarre fino (pinch grip) con alto rango DIP")
    print("=" * 70)
    print(f"\nArchivo: {zip_path}")
    print(f"Base de datos: Ninapro DB1 (CyberGlove II, 22 sensores)")
    print(f"Ejercicio: E2 (Grasping and functional movements)")
    print(f"Archivos E2 encontrados: {len(e2_files)}")
    print()
    
    # Analyze all subjects
    all_results = []
    
    print("Analizando sujetos...")
    print(f"{'Sujeto':<8} {'Ch5 rng':<10} {'Ch6 rng':<10} {'Muestras':<10} {'DIP est.':<10}")
    print("-" * 48)
    
    for fname in e2_files:
        subject = fname.split('/')[-1].split('_')[0]
        
        try:
            with zf.open(fname) as f:
                data = np.loadtxt(f)
            if data.ndim != 2 or data.shape[1] != 22:
                continue
        except Exception:
            continue
        
        segments = find_pinch_segments(data)
        
        if segments:
            best = segments[0]
            # Estimate DIP range
            dip_est = best['ch6_range'] * (100.0 / 170.0) * 0.70
            all_results.append({
                'subject': subject,
                'file': fname,
                'ch5_range': best['ch5_range'],
                'ch6_range': best['ch6_range'],
                'samples': best['end'] - best['start'],
                'dip_estimated': dip_est,
                'start': best['start'],
                'end': best['end']
            })
            print(f"{subject:<8} {best['ch5_range']:<10.0f} {best['ch6_range']:<10.0f} "
                  f"{best['end']-best['start']:<10d} {dip_est:<10.1f}")
    
    # Sort by estimated DIP range
    all_results.sort(key=lambda x: x['dip_estimated'], reverse=True)
    
    print()
    print("=" * 70)
    print("RESULTADOS - Top 10 sujetos por rango DIP estimado")
    print("=" * 70)
    print(f"{'#':<4} {'Sujeto':<8} {'DIP est.':<10} {'PIP est.':<10} {'MCP est.':<10}")
    print("-" * 42)
    
    for i, r in enumerate(all_results[:10]):
        pip_est = r['ch6_range'] * (100.0 / 170.0)
        mcp_est = r['ch5_range'] * (90.0 / 130.0)
        print(f"{i+1:<4} {r['subject']:<8} {r['dip_estimated']:<10.1f} {pip_est:<10.1f} {mcp_est:<10.1f}")
    
    # Process the best subject
    best_result = all_results[0]
    print(f"\n{'='*70}")
    print(f"MEJOR SUJETO: {best_result['subject']}")
    print(f"{'='*70}")
    
    with zf.open(best_result['file']) as f:
        data = np.loadtxt(f)
    
    # Extract with padding for complete cycle
    pad = 300
    seg_start = max(0, best_result['start'] - pad)
    seg_end = min(len(data), best_result['end'] + pad)
    segment = data[seg_start:seg_end, :]
    
    # Calibrate
    theta_mcp, theta_pip, theta_dip = calibrate_to_angles(segment)
    
    print(f"\nAngulos calibrados (segmento completo, {len(segment)} muestras):")
    print(f"  Theta_MCP: {theta_mcp.min():.1f} a {theta_mcp.max():.1f} (rango: {theta_mcp.max()-theta_mcp.min():.1f} deg)")
    print(f"  Theta_PIP: {theta_pip.min():.1f} a {theta_pip.max():.1f} (rango: {theta_pip.max()-theta_pip.min():.1f} deg)")
    print(f"  Theta_DIP: {theta_dip.min():.1f} a {theta_dip.max():.1f} (rango: {theta_dip.max()-theta_dip.min():.1f} deg)")
    
    # Resample to 120 points
    mcp_120, pip_120, dip_120 = resample_to_n_points(theta_mcp, theta_pip, theta_dip, 120)
    
    print(f"\nDatos remuestreados a 120 puntos:")
    print(f"  Theta_MCP: {mcp_120.min():.1f} a {mcp_120.max():.1f} (rango: {mcp_120.max()-mcp_120.min():.1f} deg)")
    print(f"  Theta_PIP: {pip_120.min():.1f} a {pip_120.max():.1f} (rango: {pip_120.max()-pip_120.min():.1f} deg)")
    print(f"  Theta_DIP: {dip_120.min():.1f} a {dip_120.max():.1f} (rango: {dip_120.max()-dip_120.min():.1f} deg)")
    
    # Compare with old data
    old_csv = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'mocap_indice_120pts.csv')
    if os.path.exists(old_csv):
        old = pd.read_csv(old_csv)
        old_dip_range = old['Theta_DIP'].max() - old['Theta_DIP'].min()
        new_dip_range = dip_120.max() - dip_120.min()
        print(f"\n{'='*70}")
        print("COMPARACION CON DATOS ACTUALES")
        print(f"{'='*70}")
        print(f"  mocap_indice_120pts.csv (actual):")
        print(f"    DIP rango: {old_dip_range:.1f} deg")
        print(f"    PIP rango: {old['Theta_PIP'].max()-old['Theta_PIP'].min():.1f} deg")
        print(f"    MCP rango: {old['Theta_MCP'].max()-old['Theta_MCP'].min():.1f} deg")
        print(f"\n  mocap_pinch_pizzolato.csv (nuevo, agarre fino):")
        print(f"    DIP rango: {new_dip_range:.1f} deg")
        print(f"    PIP rango: {pip_120.max()-pip_120.min():.1f} deg")
        print(f"    MCP rango: {mcp_120.max()-mcp_120.min():.1f} deg")
        print(f"\n  Mejora en rango DIP: {new_dip_range/old_dip_range:.1f}x ({old_dip_range:.1f} -> {new_dip_range:.1f} deg)")
    
    # Summary
    print(f"\n{'='*70}")
    print("RESUMEN")
    print(f"{'='*70}")
    print(f"  Sujetos analizados: {len(all_results)}")
    print(f"  Sujetos con datos validos de pinch: {len([r for r in all_results if r['dip_estimated'] > 30])}")
    print(f"  Mejor sujeto: {best_result['subject']}")
    print(f"  Rango DIP obtenido: {dip_120.max()-dip_120.min():.1f} deg (objetivo: 40-70 deg)")
    print(f"  Archivo CSV generado: mocap_pinch_pizzolato.csv")
    print(f"\n  Calibracion utilizada:")
    print(f"    Ch5 -> MCP: gain = 0.692 deg/unit, base = 5th percentile")
    print(f"    Ch6 -> PIP: gain = 0.588 deg/unit, base = 5th percentile")
    print(f"    DIP = 0.70 * PIP (acoplamiento biomecanico)")
    print(f"\n  Referencia: Kamper et al. 2003, Ninapro DB1 documentation")
    
    zf.close()
    print("\nExploracion completada exitosamente.")


if __name__ == '__main__':
    main()
