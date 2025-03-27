# Wireless Communication Project (OOK & OFDM)

## Overview

This project was developed as part of the *Wireless Communication* lab at TU Berlin. The goal was to implement and analyze two digital transmission systems using **On-Off-Keying (OOK)** and **Orthogonal Frequency Division Multiplexing (OFDM)** in **Matlab**. Both transmitter and receiver architectures were designed, simulated, and evaluated.

## Structure

- `OOK/`: Matlab scripts for OOK transmission  
- `OFDM/`: Matlab scripts for OFDM transmission  

---

## OOK (On-Off-Keying)

A basic digital transmission system using **2-ASK modulation** (Amplitude Shift Keying). The data stream is scrambled, combined with a preamble, pulse-shaped, and transmitted with a DC offset.

### Key Features

- Incoherent demodulation (no frequency or phase synchronization needed)
- Upconversion to an intermediate frequency (IF)
- Optional matched filtering for improved SNR

### Transmitter

- Data source & scrambler  
- 2-ASK modulation  
- FIR pulse shaping filter  
- DC offset addition  
- Upconversion to IF  

### Receiver

- Downconversion to baseband  
- Magnitude detection  
- DC offset removal  
- Symbol timing synchronization via preamble correlation  
- Optional matched filter  
- Decision logic & descrambler  

---

## OFDM (Orthogonal Frequency Division Multiplexing)

A more advanced system based on the IEEE 802.11a WLAN standard, using 64 subcarriers to transmit data in parallel.

### Key Features

- Supports BPSK, QPSK, 16QAM, and 64QAM  
- IFFT/FFT for modulation and demodulation  
- Guard intervals to reduce inter-symbol interference (ISI)  
- Channel and phase correction using long preamble and pilot tones  

### Transmitter

- SYNC generator (short and long preambles)  
- SIGNAL field generation  
- Modulation mapping  
- Pilot insertion  
- IFFT and guard interval insertion  
- Frame assembly and transmission  

### Receiver

- Frame detection using Schmidl-Cox correlation  
- Frequency offset correction  
- Guard interval removal  
- Channel estimation and equalization  
- Phase correction  
- Demapping & descrambling  

---

## Authors

- Anton Valentin Dilg  
- Niclas Eric Schenk  
- Ramon Rennert  
- Ramin Leon Neymeyer  
