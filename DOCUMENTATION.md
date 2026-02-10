# Dokumentacja projektu dna_encoder

## Spis treści

1. [Przegląd projektu](#przegląd-projektu)
2. [Architektura systemu](#architektura-systemu)
3. [Schemat kodowania bitów](#schemat-kodowania-bitów)
4. [Algorytm backtrackingu](#algorytm-backtrackingu)
5. [System walidacji biochemicznej](#system-walidacji-biochemicznej)
6. [Konfiguracja i profile](#konfiguracja-i-profile)
7. [Interfejs użytkownika](#interfejs-użytkownika)
8. [Szczegóły implementacyjne](#szczegóły-implementacyjne)
9. [Przykłady użycia](#przykłady-użycia)
10. [Zastosowania w bezpieczeństwie DNA](#zastosowania-w-bezpieczeństwie-dna)

---

## Przegląd projektu

### Cel projektu

Projekt `dna_encoder` to system kodowania danych binarnych i tekstowych w sekwencjach DNA z kontrolą jakości biochemicznej. W odróżnieniu od `dna_generator`, który generuje losowe sekwencje DNA spełniające kryteria jakości, `dna_encoder` umożliwia **deterministyczne kodowanie konkretnych danych** w sekwencjach DNA, które można następnie zdekodować z powrotem do oryginalnej formy.

### Główne funkcjonalności

- **Kodowanie tekstu do DNA**: Konwersja dowolnego tekstu UTF-8 na sekwencje DNA
- **Kodowanie binarne**: Bezpośrednie kodowanie ciągów bitów
- **Dekodowanie**: Odzyskiwanie oryginalnych danych z sekwencji DNA
- **Kontrola jakości biochemicznej**: Walidacja zawartości GC, temperatury topnienia, struktur wtórnych
- **Algorytm backtrackingu**: Inteligentny wybór par nukleotydów spełniających ograniczenia
- **Profile walidacji**: Konfigurowalne poziomy rygorystyczności (strict, pcr_friendly, relaxed, sequence_only)
- **Wsparcie dla primerów PCR**: 500+ zwalidowanych primerów w zestawie

### Zastosowania

- **Przechowywanie danych w DNA**: Archiwizacja danych cyfrowych w syntetycznym DNA
- **Znakowanie biologiczne**: Ukryte znaczniki w organizmach genetycznie modyfikowanych
- **Bezpieczeństwo tagów DNA**: Analiza i weryfikacja autentyczności sekwencji
- **Steganografia molekularna**: Ukrywanie informacji w sekwencjach DNA
- **Badania nad kodowaniem DNA**: Eksperymenty z różnymi schematami kodowania

---

## Architektura systemu

### Komponenty główne

```
dna_encoder/
├── encoder.py           # Główna klasa DNAEncoder
├── decoder.py           # Klasa DNADecoder do dekodowania
├── backtracking_engine.py  # Algorytm backtrackingu
├── config.py            # Konfiguracja DNAEncoderConfig
├── profile_loader.py    # Ładowanie profili z JSON
├── results.py           # Struktura wyników DNAResult
├── utils.py             # Funkcje pomocnicze
├── exceptions.py        # Hierarchia wyjątków
├── default_profiles.json   # Domyślne profile walidacji
└── __main__.py          # Interfejs CLI
```

### Zależności zewnętrzne

```
dna_commons>=0.1.0       # Wspólna biblioteka DNA
├── DNAValidator         # Walidacja sekwencji
├── Primer3Adapter       # Obliczenia termodynamiczne
├── DeterministicRandom  # Kontrolowany RNG
└── SequenceAnalyzer     # Analiza sekwencji
```

### Separacja odpowiedzialności

#### 1. **Warstwa API** (`encoder.py`, `decoder.py`)
- `DNAEncoder`: Główna klasa kodowania
- `DNADecoder`: Klasa dekodowania
- Obsługa błędów wysokiego poziomu
- Logowanie i metryki

#### 2. **Warstwa algorytmiczna** (`backtracking_engine.py`)
- `BacktrackingEngine`: Implementacja algorytmu
- Heurystyki wyboru par nukleotydów
- Statystyki backtrackingu

#### 3. **Warstwa konfiguracji** (`config.py`, `profile_loader.py`)
- `DNAEncoderConfig`: Centralna konfiguracja
- `ProfileLoader`: Ładowanie profili z JSON
- Walidacja parametrów

#### 4. **Warstwa walidacji** (z `dna_commons`)
- `DNAValidator`: Walidacja kryteriów biochemicznych
- `Primer3Adapter`: Obliczenia termodynamiczne

### Przepływ danych - Kodowanie

```
Tekst/Bity → string_to_bits() → BacktrackingEngine → DNAValidator → Sekwencja DNA
                    ↓                    ↓                  ↓              ↓
            "01001000..."    → wybór par nukleotydów → walidacja → "ATGC..."
```

### Przepływ danych - Dekodowanie

```
Sekwencja DNA → _extract_encoded_part() → _decode_pairs_to_bits() → bits_to_string() → Tekst
    ↓                    ↓                        ↓                      ↓
"ATGC..."    →    "ATAT..."       →          "01001000..."    →     "Hello"
```

---

## Schemat kodowania bitów

### Podstawy teoretyczne

System używa **16 unikalnych par nukleotydów** do reprezentacji bitów. Każdy bit (`0` lub `1`) może być zakodowany za pomocą jednej z 8 dostępnych par nukleotydów:

```python
BIT_TO_NUCLEOTIDES = {
    '0': ['TA', 'TT', 'GC', 'CC', 'AC', 'AG', 'GT', 'CT'],
    '1': ['AT', 'AA', 'CG', 'GG', 'TC', 'TG', 'GA', 'CA']
}
```

### Właściwości schematu

#### 1. **Redundancja kodowania**
Każdy bit ma 8 alternatywnych reprezentacji, co umożliwia:
- Wybór pary najlepiej pasującej do ograniczeń biochemicznych
- Optymalizację zawartości GC
- Unikanie struktur wtórnych

#### 2. **Jednoznaczność dekodowania**
Każda para nukleotydów jednoznacznie określa bit:
```
TA, TT, GC, CC, AC, AG, GT, CT → 0
AT, AA, CG, GG, TC, TG, GA, CA → 1
```

#### 3. **Zbalansowane mapowanie**
- 8 par dla bitu `0`, 8 par dla bitu `1`
- Równomierna dystrybucja nukleotydów A, T, G, C

### Konwersja tekstu na bity

```python
def string_to_bits(text: str, encoding: str = 'utf-8') -> str:
    """Konwertuje tekst na reprezentację binarną (8 bitów na znak)."""
    return ''.join(format(byte, '08b') for byte in text.encode(encoding))

# Przykład:
# "H" → 01001000
# "Hi" → 0100100001101001
```

### Struktura zakodowanej sekwencji

```
[Primer/Initial Sequence][Zakodowane dane]
        ↓                        ↓
   "ATGCATGC..."          2 nukleotydy/bit
   (20-30 nt)             (długość = bity × 2)
```

**Przykład:**
```
Tekst: "Hi" (2 znaki)
Bity: 16 bitów (2 × 8)
Pary nukleotydów: 16 par
Długość zakodowanych danych: 32 nukleotydy
Całkowita sekwencja: primer (20 nt) + dane (32 nt) = 52 nt
```

---

## Algorytm backtrackingu

### Podstawy teoretyczne

Algorytm backtrackingu przeszukuje przestrzeń możliwych sekwencji DNA, wybierając optymalne pary nukleotydów dla każdego bitu. Gdy wybrana para prowadzi do naruszenia ograniczeń biochemicznych, algorytm cofa się i próbuje alternatywnej pary.

### Struktura przestrzeni poszukiwań

Dla każdego bitu do zakodowania:
```
Bit 0:           Bit 1:           ...
   |                |
   ├─ TA            ├─ AT
   ├─ TT            ├─ AA
   ├─ GC            ├─ CG
   ├─ CC            ├─ GG
   ├─ AC            ├─ TC
   ├─ AG            ├─ TG
   ├─ GT            ├─ GA
   └─ CT            └─ CA
```

### Implementacja techniczna

#### Główna pętla algorytmu

```python
def _execute_encoding_algorithm(self, initial_sequence, bit_sequence, base_seed, rng, stats):
    nucleotide_pairs = [BIT_TO_NUCLEOTIDES[bit].copy() for bit in bit_sequence]
    sequence_parts = [initial_sequence]
    available_options = []

    if nucleotide_pairs:
        available_options.append(nucleotide_pairs[0].copy())

    position = 0
    backtrack_count = 0

    while position < len(bit_sequence):
        stats['total_attempts'] += 1

        # Sprawdź czy są dostępne opcje
        if not available_options or not available_options[-1]:
            backtrack_result = self._handle_backtrack(...)
            if backtrack_result is None:
                return ''.join(sequence_parts), backtrack_count
            position, backtrack_count = backtrack_result
            continue

        # Wybierz i waliduj parę
        chosen_pair = self._choose_and_validate_pair(...)

        if chosen_pair:
            sequence_parts.append(chosen_pair)
            position += 1
            if position < len(bit_sequence):
                available_options.append(nucleotide_pairs[position].copy())

    return ''.join(sequence_parts), backtrack_count
```

#### Mechanizm cofania

```python
def _handle_backtrack(self, sequence_parts, available_options, position,
                      backtrack_count, max_backtracks, stats):
    # Nie cofaj się poza sekwencję początkową
    if len(sequence_parts) <= 1:
        return None  # Niemożliwe zakodowanie

    # Cofnij ostatnią parę
    sequence_parts.pop()
    available_options.pop()
    position -= 1
    backtrack_count += 1

    # Sprawdź limit backtrackingu
    if backtrack_count > max_backtracks:
        return None

    return position, backtrack_count
```

### System heurystyk

#### Punktacja kandydatów

System używa wielokryterialnej funkcji oceny do wyboru optymalnej pary nukleotydów:

```python
def _calculate_heuristic_score(self, pair, window, target_gc, recent_context,
                               candidate_sequence, last_pair, last_pair_complement,
                               gc_weight, diversity_weight, pair_repeat_weight,
                               pair_complement_weight, novelty_weight):

    # 1. Odległość od docelowego GC (tylko jeśli walidacja włączona)
    gc_score = 0.0
    if target_gc is not None:
        gc_fraction = (window.count('G') + window.count('C')) / len(window)
        gc_score = gc_weight * abs(gc_fraction - target_gc)

    # 2. Kara za brak różnorodności nukleotydów
    diversity_score = (window.count(pair[0]) + window.count(pair[1])) / len(window)
    diversity_contribution = diversity_weight * diversity_score

    # 3. Kara za powtórzenie pary
    pair_repeat_score = window.count(pair) / max(1, len(window) // 2)
    immediate_repeat = 1.0 if last_pair and pair == last_pair else 0.0
    repeat_contribution = pair_repeat_weight * (pair_repeat_score + immediate_repeat)

    # 4. Kara za komplementarność
    complement_contribution = pair_complement_weight if pair == last_pair_complement else 0.0

    # 5. Kara za okresowość (anty-nowość)
    novelty_hits = 0
    for k in (4, 6, 8):
        if len(candidate_sequence) >= k:
            kmer = candidate_sequence[-k:]
            novelty_hits += recent_context.count(kmer)
    novelty_contribution = novelty_weight * novelty_hits

    return gc_score + diversity_contribution + repeat_contribution + \
           complement_contribution + novelty_contribution
```

#### Kary twarde (Hard penalties)

Stosowane gdy kandydat narusza twarde ograniczenia:

```python
def _calculate_hard_penalties(self, window):
    hard_penalty = 0.0

    if self.validator.rules.homopolymer_runs:
        if not self.validator._check_homopolymer_runs(window)[0]:
            hard_penalty += 300_000.0  # Ciągi jednakowych nukleotydów

    if self.validator.rules.three_prime_stability:
        if not self.validator._check_3_prime_stability(window):
            hard_penalty += 200_000.0  # Niestabilność końca 3'

    if self.validator.rules.dinucleotide_repeats:
        if not self.validator._check_dinucleotide_repeats(window)[0]:
            hard_penalty += 150_000.0  # Powtórzenia dinukleotydów

    return hard_penalty
```

#### Domyślne wagi heurystyk

```python
heuristic_weights = {
    'gc_balance': 1.0,       # Waga dla równowagi GC
    'diversity': 0.05,       # Waga dla różnorodności nukleotydów
    'pair_repeat': 0.3,      # Waga dla unikania powtórzeń par
    'pair_complement': 0.2,  # Waga dla unikania komplementów
    'novelty': 0.01          # Waga dla nowości (anty-okresowość)
}
```

### Tryby wyboru

#### Tryb deterministyczny
```python
if self.config.enable_deterministic:
    # Wybierz deterministycznie spośród najlepszych
    tolerance = 1e-6
    best_candidates = [pair for score, pair in scored
                       if abs(score - best_score) <= tolerance]
    return self._select_with_deterministic_tiebreak(best_candidates, base_seed, position, rng)
```

#### Tryb losowy
```python
if not self.config.enable_deterministic:
    # Wybierz losowo z marginesem
    margin = 0.02
    eligible = [pair for score, pair in scored if score <= best_score + margin]
    return rng.choice(eligible)
```

---

## System walidacji biochemicznej

### Metryki jakości DNA

System walidacji sprawdza następujące kryteria biochemiczne (dziedziczone z `dna_commons`):

#### 1. **Zawartość GC**
- **Zakres**: Kontrolowany parametrami `min_gc` i `max_gc`
- **Domyślnie**: 45-55% (strict) lub 40-60% (relaxed)
- **Znaczenie**: Wpływa na stabilność termiczną duplex DNA

```python
gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
```

#### 2. **Temperatura topnienia (Tm)**
- **Zakres**: Kontrolowany parametrami `min_tm` i `max_tm`
- **Domyślnie**: 55-65°C
- **Metoda**: Algorytm nearest-neighbor (primer3) lub fallback Wallace

#### 3. **Struktury hairpin (spinki do włosów)**
```
5'-ATGCGCATGCGCAT-3'
    |||||||
    TACGCGTAT  ← pętla hairpin
```
- **Próg**: `max_hairpin_tm` (domyślnie 25-32°C)
- **Problem**: Może blokować amplifikację PCR

#### 4. **Homodimery (dupleksy)**
```
5'-ATGCGCAT-3'
   ||||||||
3'-TACGCGTA-5'  ← homodimer
```
- **Próg**: `max_homodimer_tm` (domyślnie 25-32°C)
- **Problem**: Self-priming, niższa efektywność PCR

#### 5. **Homopolymer runs (ciągi jednakowych nukleotydów)**
```
ATGCCCCCATGC  ← 5× C (za długie)
ATGCCCATGC    ← 3× C (akceptowalne)
```
- **Limit**: `max_homopolymer_length` (domyślnie 3-4)
- **Problem**: Błędy polimerazy, problemy z sekwencjonowaniem

#### 6. **Powtórzenia dinukleotydowe**
```
ATGCATATATATGC  ← 4× AT (za dużo)
ATGCATATGC      ← 2× AT (akceptowalne)
```
- **Limit**: `max_dinucleotide_repeats` (domyślnie 2-3)
- **Problem**: Niestabilność strukturalna

#### 7. **Stabilność końca 3'**
```
...ATGCGCGC-3'  ← 4× GC na końcu (za stabilne)
...ATGCATGC-3'  ← 3× GC (akceptowalne)
```
- **Limit**: `max_3prime_gc` (domyślnie 3-4 nukleotydy GC w ostatnich 5 nt)
- **Problem**: Niespecyficzne priming w PCR

### Okno analizy

Walidacja używa przesuwnego okna do analizy lokalnych właściwości:

```python
def _validate_candidate(self, candidate_sequence, stats):
    window_size = self.config.window_size  # Domyślnie 20

    if len(candidate_sequence) >= window_size:
        window_seq = candidate_sequence[-window_size:]
    else:
        window_seq = candidate_sequence

    return self.validator.validate_window(window_seq)
```

---

## Konfiguracja i profile

### DNAEncoderConfig

Centralna klasa konfiguracji systemu:

```python
@dataclass
class DNAEncoderConfig:
    # Parametry primer3 (termodynamika)
    primer3_mv_conc: float = 50.0      # Mg2+ [mM]
    primer3_dv_conc: float = 4.0       # Divalent cations [mM]
    primer3_dntp_conc: float = 0.5     # dNTP [mM]
    primer3_dna_conc: float = 50.0     # DNA [nM]
    primer3_temp_c: float = 37.0       # Temperatura [°C]

    # Ograniczenia jakości sekwencji
    min_gc: float = 0.45
    max_gc: float = 0.55
    min_tm: float = 55.0
    max_tm: float = 65.0
    max_hairpin_tm: float = 30.0
    max_homodimer_tm: float = 30.0

    # Parametry algorytmu
    window_size: int = 20
    max_homopolymer_length: int = 4
    max_dinucleotide_repeats: int = 2
    max_3prime_gc: int = 3
    max_backtrack_attempts: int = 1000
    enable_backtrack_heuristics: bool = True

    # Tryb deterministyczny
    enable_deterministic: bool = True
    default_seed: int = None

    # Profil walidacji
    validation_profile: str = None
```

### Profile walidacji

Profile są zdefiniowane w `default_profiles.json`:

#### 1. **Profile `strict`**
```json
{
  "description": "Strict quality constraints for high-quality sequences",
  "rules": {
    "gc_content": true,
    "melting_temperature": true,
    "hairpin_structures": true,
    "homodimer_structures": true,
    "homopolymer_runs": true,
    "dinucleotide_repeats": true,
    "three_prime_stability": true
  },
  "params": {
    "min_gc": 0.45, "max_gc": 0.55,
    "min_tm": 57.0, "max_tm": 63.0,
    "max_hairpin_tm": 25.0,
    "max_homodimer_tm": 25.0,
    "max_homopolymer_length": 3,
    "max_dinucleotide_repeats": 2,
    "max_3prime_gc": 3
  }
}
```

#### 2. **Profile `pcr_friendly`**
```json
{
  "description": "Profile optimized for generating long PCR-compatible sequences",
  "rules": { /* wszystkie reguły włączone */ },
  "params": {
    "min_gc": 0.40, "max_gc": 0.60,
    "min_tm": 54.0, "max_tm": 66.0,
    "max_hairpin_tm": 32.0,
    "max_homodimer_tm": 32.0,
    "max_homopolymer_length": 4,
    "max_dinucleotide_repeats": 3,
    "max_3prime_gc": 4
  }
}
```

#### 3. **Profile `relaxed`**
```json
{
  "description": "Relaxed constraints for general-purpose use",
  "params": {
    "min_gc": 0.40, "max_gc": 0.60,
    "min_tm": 52.0, "max_tm": 68.0,
    /* szersze zakresy */
  }
}
```

#### 4. **Profile `sequence_only`**
```json
{
  "description": "Only basic sequence quality checks (no thermodynamic validation)",
  "rules": {
    "gc_content": false,
    "melting_temperature": false,
    "hairpin_structures": false,
    "homodimer_structures": false,
    "homopolymer_runs": true,
    "dinucleotide_repeats": true,
    "three_prime_stability": true
  }
}
```

#### 5. **Profile `none`**
```json
{
  "description": "No validation rules - completely random generation",
  "rules": { /* wszystkie false */ },
  "params": { /* brak ograniczeń */ }
}
```

### Hierarchia konfiguracji

1. **Wartości domyślne** w `DNAEncoderConfig`
2. **Profile walidacji** z plików JSON
3. **Parametry użytkownika** (przesłaniają profil)

```python
# Kolejność ładowania:
config = DNAEncoderConfig(
    validation_profile='strict',  # Najpierw ładowany profil
    min_gc=0.48                   # Potem nadpisanie parametru
)
```

### Własne profile użytkownika

Użytkownicy mogą definiować własne profile w pliku `user_profiles.json`:

```json
{
  "profiles": {
    "my_custom_profile": {
      "description": "My custom validation settings",
      "rules": {
        "gc_content": true,
        "melting_temperature": false,
        /* ... */
      },
      "params": {
        "min_gc": 0.42,
        "max_gc": 0.58,
        /* ... */
      }
    }
  }
}
```

---

## Interfejs użytkownika

### API programistyczne

#### Podstawowe użycie

```python
from dna_encoder import DNAEncoder, DNADecoder

# Inicjalizacja z domyślną konfiguracją
encoder = DNAEncoder()
decoder = DNADecoder()

# Primer (sekwencja początkowa)
primer = "CATCTATCCCTTCGAACGAC"

# Kodowanie tekstu
result = encoder.encode_text(primer, "Hello, DNA!")

if result.success:
    print(f"Sekwencja: {result.sequence}")
    print(f"Długość: {len(result.sequence)} nt")

    # Dekodowanie
    decoded = decoder.decode_sequence(result.sequence, primer)
    print(f"Zdekodowano: {decoded}")
else:
    print(f"Błąd: {result.error_message}")
```

#### Kodowanie bitów

```python
# Bezpośrednie kodowanie ciągu bitów
bits = "01001000"  # = 'H' w ASCII
result = encoder.encode_bits(primer, bits)

if result.success:
    decoded = decoder.decode_sequence(result.sequence, primer)
    print(f"Zdekodowano: {decoded}")  # "H"
```

#### Zaawansowana konfiguracja

```python
from dna_encoder import DNAEncoderConfig, DNAEncoder

# Konfiguracja z profilem
config = DNAEncoderConfig(
    validation_profile='strict',
    max_backtrack_attempts=2000,
    window_size=25
)

encoder = DNAEncoder(config)
result = encoder.encode_text(primer, "Secret message", max_attempts=10)
```

#### Szybkie funkcje pomocnicze

```python
from dna_encoder import quick_encode, quick_decode

# Szybkie kodowanie z domyślnymi ustawieniami
result = quick_encode("Hello World!")
if result.success:
    text = quick_decode(result.sequence)
```

### Interfejs CLI

#### Podstawowe polecenia

```bash
# Kodowanie tekstu
python -m dna_encoder --encode "Hello" --initial ATGCATGC

# Kodowanie bitów
python -m dna_encoder --encode-bits "01001100" --initial ATGC

# Dekodowanie sekwencji
python -m dna_encoder --decode "ATGCATGCAAGGTTCC..." --initial ATGCATGC

# Wyświetlenie informacji
python -m dna_encoder

# Uruchomienie testów
python -m dna_encoder --test
```

#### Pełna lista opcji CLI

```bash
Argumenty główne:
  --encode TEXT         Zakoduj podany tekst do DNA
  --encode-bits BITS    Zakoduj ciąg bitów do DNA
  --decode DNA          Zdekoduj sekwencję DNA
  --initial SEQ         Sekwencja początkowa (primer)

Opcje trybu:
  --max-attempts N      Limit prób kodowania (domyślnie 5)
  --random              Wyłącz tryb deterministyczny
  --seed N              Wymuś konkretny seed
  --no-heuristics       Wyłącz heurystyki wyboru par

Profile walidacji:
  --profile NAME        Wybierz profil (strict, pcr_friendly, relaxed, sequence_only, none)
  --profile-file JSON   Ścieżka do własnego profilu JSON

Parametry biochemiczne:
  --min-gc, --max-gc    Zakres zawartości GC (0.0-1.0)
  --min-tm, --max-tm    Zakres temperatury topnienia (°C)
  --max-hairpin         Limit Tm hairpinów (°C)
  --max-homodimer       Limit Tm homodimerów (°C)
  --window-size         Rozmiar okna analizy

Wyłączanie reguł:
  --no-gc-check         Wyłącz kontrolę GC
  --no-tm-check         Wyłącz kontrolę Tm
  --no-hairpin-check    Wyłącz kontrolę hairpinów
  --no-homodimer-check  Wyłącz kontrolę homodimerów
  --no-homopolymer-check  Wyłącz kontrolę homopolimerów
  --no-dinucleotide-check Wyłącz kontrolę dinukleotydów
  --no-3prime-check     Wyłącz kontrolę końca 3'

Wyjście:
  --output FILE         Zapisz wynik do pliku
  --format FMT          Format wyjścia (text, json, fasta)
  --sequences-only      Tylko sekwencje (bez statystyk)
  --csv-file FILE       Eksportuj analizę sliding window do CSV
  --quiet               Ogranicz logowanie
  --verbose             Szczegółowe logowanie
```

#### Przykłady wyjścia

**Format text (domyślny):**
```
Encoding text: 'Hello'
Initial sequence: ATGCATGC
✓ Encoding completed successfully
     1: ATGCATGC ATCGTAGC TAGCTAGC ATGC
Statistics:
  Attempts: 1
  Backtracks: 12
  Total choices: 47
  Max depth: 52
```

**Format JSON:**
```json
{
  "success": true,
  "sequence": "ATGCATGCATCGTAGCTAGCTAGCATGC...",
  "encoding_stats": {
    "attempts": 1,
    "backtrack_count": 12,
    "total_attempts": 47,
    "max_depth_reached": 52
  }
}
```

**Format FASTA:**
```
>encoded|success=1
ATGCATGCATCGTAGCTAGCTAGCATGC...
```

---

## Szczegóły implementacyjne

### Struktura wyników DNAResult

```python
@dataclass
class DNAResult:
    success: bool                           # Czy kodowanie się powiodło
    sequence: Optional[str] = None          # Wygenerowana sekwencja DNA
    error_message: Optional[str] = None     # Komunikat błędu
    encoding_stats: Optional[Dict] = None   # Statystyki kodowania
    validation_details: Optional[Dict] = None  # Szczegóły walidacji

    def get_sequence_length(self) -> int: ...
    def get_encoding_attempts(self) -> int: ...
    def get_backtrack_count(self) -> int: ...
```

### Statystyki kodowania

```python
encoding_stats = {
    'seed': 42,                      # Użyty seed RNG
    'heuristics_enabled': True,      # Czy heurystyki włączone
    'initial_sequence_length': 20,   # Długość primera
    'target_bits': 88,               # Liczba bitów do zakodowania
    'backtrack_count': 15,           # Liczba cofnięć
    'total_attempts': 103,           # Próby dodania par
    'max_depth_reached': 108,        # Maksymalna długość sekwencji
    'validation_failures': {
        'gc_content': 3,
        'melting_temp': 2,
        'homopolymers': 5,
        'dinucleotide_repeats': 1,
        'three_prime_stability': 4,
        'hairpin': 0,
        'homodimer': 0,
    },
    'sequence_analysis': {
        'gc_content': 0.52,
        'complexity_k2': 0.94,
        'longest_homopolymer': {'nucleotide': 'A', 'length': 3}
    }
}
```

### Hierarchia wyjątków

```python
class DNAEncodingError(Exception):
    """Bazowy wyjątek dla wszystkich błędów"""
    pass

class ConfigurationError(DNAEncodingError):
    """Błędy konfiguracji"""
    pass

class ValidationError(DNAEncodingError):
    """Błędy walidacji sekwencji"""
    def __init__(self, message, sequence=None, failed_checks=None): ...

class BacktrackingError(DNAEncodingError):
    """Algorytm nie znalazł rozwiązania"""
    def __init__(self, message, attempts=None, last_position=None): ...

class EncodingError(DNAEncodingError):
    """Błędy w procesie kodowania"""
    def __init__(self, message, input_data=None, bit_position=None): ...

class DecodingError(DNAEncodingError):
    """Błędy w procesie dekodowania"""
    def __init__(self, message, dna_sequence=None, invalid_pair=None): ...

class InputError(DNAEncodingError):
    """Błędy danych wejściowych"""
    pass

class Primer3Error(DNAEncodingError):
    """Błędy biblioteki primer3"""
    pass
```

### Deterministyczny RNG

System używa deterministycznego generatora liczb pseudolosowych dla powtarzalności:

```python
class DeterministicRandom:
    def __init__(self, seed=None, deterministic=True):
        self._random = random.Random()
        if deterministic:
            self._random.seed(seed or 42)
```

### Generowanie seeda z parametrów

```python
def generate_seed_from_string(initial_sequence: str, additional_data: str = "") -> int:
    """Generuje deterministyczny seed z parametrów wejściowych."""
    combined = f"{initial_sequence}:{additional_data}"
    hash_object = hashlib.md5(combined.encode())
    return int(hash_object.hexdigest()[:8], 16)
```

---

## Przykłady użycia

### 1. Podstawowe kodowanie i dekodowanie

```python
from dna_encoder import DNAEncoder, DNADecoder

encoder = DNAEncoder()
decoder = DNADecoder()

primer = "CATCTATCCCTTCGAACGAC"
message = "Hello, DNA!"

# Kodowanie
result = encoder.encode_text(primer, message)

if result.success:
    print(f"Zakodowano: {result.sequence}")
    print(f"Długość: {len(result.sequence)} nt")

    # Dekodowanie
    decoded = decoder.decode_sequence(result.sequence, primer)
    print(f"Zdekodowano: {decoded}")
    assert decoded == message
```

### 2. Użycie profili walidacji

```python
from dna_encoder import DNAEncoder, DNAEncoderConfig

# Różne profile dla różnych zastosowań
profiles = ['sequence_only', 'relaxed', 'pcr_friendly', 'strict']

primer = "CATCTATCCCTTCGAACGAC"
message = "Test"

for profile in profiles:
    config = DNAEncoderConfig(validation_profile=profile)
    encoder = DNAEncoder(config)
    result = encoder.encode_text(primer, message)

    status = "✓" if result.success else "✗"
    backtracks = result.encoding_stats.get('backtrack_count', 0) if result.success else "N/A"
    print(f"{profile:15} {status} backtracks: {backtracks}")
```

### 3. Kodowanie danych binarnych

```python
from dna_encoder import DNAEncoder, DNADecoder

encoder = DNAEncoder()
decoder = DNADecoder()

primer = "ATGCATGCATGCATGCATGC"

# Kodowanie bajtów (np. obrazu, pliku)
binary_data = bytes([0x48, 0x65, 0x6C, 0x6C, 0x6F])  # "Hello"
bits = ''.join(format(byte, '08b') for byte in binary_data)

result = encoder.encode_bits(primer, bits)

if result.success:
    # Dekodowanie
    decoded_text = decoder.decode_sequence(result.sequence, primer)
    decoded_bytes = decoded_text.encode('utf-8')
    print(f"Zdekodowano: {decoded_bytes}")  # b'Hello'
```

### 4. Użycie zwalidowanych primerów

```python
from pathlib import Path
from dna_encoder import DNAEncoder, DNAEncoderConfig

# Wczytaj primery z pliku
primers_file = Path('primers.txt')
with open(primers_file) as f:
    primers = [line.strip() for line in f if line.strip()]

print(f"Załadowano {len(primers)} primerów")

# Użyj losowego primera
import random
primer = random.choice(primers)

config = DNAEncoderConfig(validation_profile='pcr_friendly')
encoder = DNAEncoder(config)

result = encoder.encode_text(primer, "DNA Storage Test")
if result.success:
    print(f"Primer: {primer}")
    print(f"Sekwencja: {result.sequence}")
```

### 5. Analiza statystyk kodowania

```python
import statistics
from dna_encoder import DNAEncoder, DNAEncoderConfig

config = DNAEncoderConfig(validation_profile='pcr_friendly')
encoder = DNAEncoder(config)

primer = "CATCTATCCCTTCGAACGAC"
messages = ["Hi", "Test", "DNA", "Hello World"]

all_backtracks = []
all_times = []

import time
for msg in messages:
    start = time.time()
    result = encoder.encode_text(primer, msg, max_attempts=10)
    elapsed = time.time() - start

    if result.success:
        backtracks = result.encoding_stats.get('backtrack_count', 0)
        all_backtracks.append(backtracks)
        all_times.append(elapsed)
        print(f"'{msg}': {backtracks} backtracks, {elapsed*1000:.1f}ms")

print(f"\nStatystyki:")
print(f"  Średnie backtracking: {statistics.mean(all_backtracks):.1f}")
print(f"  Średni czas: {statistics.mean(all_times)*1000:.1f}ms")
```

### 6. Eksport analizy do CSV

```bash
# Z linii poleceń
python -m dna_encoder --encode "Test data" --initial ATGCATGC \
    --profile pcr_friendly --csv-file analysis.csv
```

```python
# Programowo
import csv
from dna_encoder import DNAEncoder, DNAEncoderConfig

config = DNAEncoderConfig(validation_profile='pcr_friendly')
encoder = DNAEncoder(config)

primer = "CATCTATCCCTTCGAACGAC"
result = encoder.encode_text(primer, "Test")

if result.success:
    sequence = result.sequence
    window_size = config.window_size

    # Analiza sliding window
    rows = []
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i+window_size]
        metrics = encoder.validator.validate_sequence(window)
        rows.append({
            'position': i,
            'window': window,
            'gc_content': metrics.gc_content,
            'tm': metrics.melting_temperature,
            'is_valid': metrics.is_valid
        })

    # Zapis do CSV
    with open('analysis.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['position', 'window', 'gc_content', 'tm', 'is_valid'])
        writer.writeheader()
        writer.writerows(rows)
```

---

## Zastosowania w bezpieczeństwie DNA

### Analiza tagów DNA

Projekt `dna_encoder` może być użyty do analizy bezpieczeństwa tagów DNA w kontekście:

#### 1. **Weryfikacja autentyczności**
```python
from dna_encoder import DNADecoder

decoder = DNADecoder()

# Znany primer producenta
authentic_primer = "CATCTATCCCTTCGAACGAC"

# Sekwencja z próbki
sample_sequence = "CATCTATCCCTTCGAACGACATAT..."

try:
    decoded = decoder.decode_sequence(sample_sequence, authentic_primer)
    print(f"Tag autentyczny: {decoded}")
except Exception as e:
    print(f"Tag nieautentyczny lub uszkodzony: {e}")
```

#### 2. **Analiza odporności na mutacje**
```python
import random

def simulate_mutation(sequence, mutation_rate=0.01):
    """Symuluje losowe mutacje w sekwencji."""
    nucleotides = list(sequence)
    for i in range(len(nucleotides)):
        if random.random() < mutation_rate:
            nucleotides[i] = random.choice('ATGC')
    return ''.join(nucleotides)

# Test odporności na mutacje
original_sequence = result.sequence
mutated_sequence = simulate_mutation(original_sequence, 0.02)

try:
    decoded = decoder.decode_sequence(mutated_sequence, primer)
    print(f"Dekodowanie po mutacjach: {decoded}")
except:
    print("Dekodowanie nieudane - zbyt wiele mutacji")
```

#### 3. **Analiza przestrzeni kodowej**
```python
# Liczba możliwych tagów dla danej długości
def calculate_tag_space(bits_length):
    """Oblicza liczbę możliwych sekwencji dla danej liczby bitów."""
    # Każdy bit ma 8 możliwych reprezentacji
    return 8 ** bits_length

# Dla 128-bitowego tagu (16 bajtów)
print(f"Przestrzeń 128-bit: {calculate_tag_space(128):.2e} możliwych sekwencji")
```

#### 4. **Steganografia molekularna**
```python
# Ukrycie danych w "losowej" sekwencji DNA
secret_message = "SECRET"
cover_primer = "ATGCATGCATGCATGCATGC"  # Wygląda jak zwykły primer

config = DNAEncoderConfig(
    validation_profile='pcr_friendly',
    enable_deterministic=True,
    default_seed=42  # Tajny klucz
)
encoder = DNAEncoder(config)

result = encoder.encode_text(cover_primer, secret_message)
# Sekwencja wygląda jak typowa sekwencja DNA
# Tylko znający primer i seed może ją zdekodować
```

---

## Załączniki

### A. Porównanie z dna_generator

| Aspekt | dna_generator | dna_encoder |
|--------|--------------|-------------|
| **Cel** | Generowanie sekwencji DNA | Kodowanie danych w DNA |
| **Wejście** | Długość sekwencji | Tekst/bity + primer |
| **Wyjście** | Losowa sekwencja | Deterministyczna sekwencja |
| **Odwracalność** | Nie | Tak (dekodowanie) |
| **Algorytm** | Backtracking (wybór 1 z 4) | Backtracking (wybór 1 z 8) |
| **Zastosowanie** | Projektowanie primerów | Przechowywanie danych |

### B. Złożoność obliczeniowa

- **Najgorszy przypadek**: O(8^n) gdzie n = liczba bitów
- **Średni przypadek**: O(n × k) gdzie k = średnia liczba backtrack'ów
- **Pamięć**: O(n) dla stosu stanów

### C. Wydajność typowa

```
Długość tekstu | Bity | Czas [ms] | Backtracks | Sukces (pcr_friendly)
1 znak        |   8  |    5-20   |     2-10   | >95%
5 znaków      |  40  |   20-100  |    10-50   | >90%
10 znaków     |  80  |  50-300   |   30-150   | >85%
20 znaków     | 160  | 100-800   |  100-500   | >75%
```

### D. Referencje

- **Primer3**: Koressaar T, Remm M. 2007. "Enhancements and modifications of primer design program Primer3"
- **DNA data storage**: Church GM et al. 2012. "Next-Generation Digital Information Storage in DNA"
- **Nearest-neighbor thermodynamics**: SantaLucia J. 1998

---

*Dokumentacja wygenerowana dla projektu dna_encoder v0.1.0*
