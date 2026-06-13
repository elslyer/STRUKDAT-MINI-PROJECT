import os
import csv
import io
import base64

from flask import Flask, render_template, request, send_file
from Bio import SeqIO
import matplotlib.pyplot as plt

app = Flask(__name__)

app.config['UPLOAD_FOLDER'] = 'uploads'

os.makedirs('uploads', exist_ok=True)
os.makedirs('results', exist_ok=True)


def calculate_gc(sequence):
    g = sequence.count('G')
    c = sequence.count('C')

    if len(sequence) == 0:
        return 0

    return ((g + c) / len(sequence)) * 100


def get_nucleotide_freq(sequence):
    return {
        'A': sequence.count('A'),
        'C': sequence.count('C'),
        'G': sequence.count('G'),
        'T': sequence.count('T')
    }


def create_plot(data):

    if len(data) == 0:
        return None

    top_plot = data[:20]

    ids = [item['id'][:12] for item in top_plot]
    gc_values = [item['gc'] for item in top_plot]

    plt.figure(figsize=(12, 6))
    plt.bar(ids, gc_values)

    plt.xlabel("Sequence ID")
    plt.ylabel("GC Content (%)")
    plt.title("Top 20 GC Content")
    plt.xticks(rotation=45)

    plt.tight_layout()

    img = io.BytesIO()

    plt.savefig(img, format='png')
    img.seek(0)

    plot_url = base64.b64encode(
        img.getvalue()
    ).decode()

    plt.close()

    return plot_url


@app.route('/', methods=['GET', 'POST'])
def index():

    if request.method == 'POST':

        files = request.files.getlist('files')

        all_sequences = []

        total_freq = {
            'A': 0,
            'C': 0,
            'G': 0,
            'T': 0
        }

        for file in files:

            if file.filename == '':
                continue

            filepath = os.path.join(
                app.config['UPLOAD_FOLDER'],
                file.filename
            )

            file.save(filepath)

            if file.filename.endswith(('.fasta', '.fa')):
                fmt = "fasta"

            elif file.filename.endswith('.fastq'):
                fmt = "fastq"

            else:
                continue

            for record in SeqIO.parse(filepath, fmt):

                seq_str = str(record.seq).upper()

                gc_content = calculate_gc(seq_str)

                freq = get_nucleotide_freq(seq_str)

                for nuc in total_freq:
                    total_freq[nuc] += freq[nuc]

                all_sequences.append({

                    'file': file.filename,
                    'id': record.id,
                    'gc': gc_content,
                    'freq': freq,
                    'length': len(seq_str)

                })

        sorted_data = sorted(
            all_sequences,
            key=lambda x: x['gc'],
            reverse=True
        )

        total_sequences = len(sorted_data)

        if total_sequences == 0:

            return render_template(
                'index.html',
                top_3=None
            )

        # ==========================
        # BASIC STATISTICS
        # ==========================

        top_3 = sorted_data[:3]
        top_10 = sorted_data[:10]

        highest_gc_seq = sorted_data[0]
        lowest_gc_seq = sorted_data[-1]

        avg_gc = sum(
            item['gc']
            for item in sorted_data
        ) / total_sequences

        avg_length = sum(
            item['length']
            for item in sorted_data
        ) / total_sequences

        max_gc = highest_gc_seq['gc']
        min_gc = lowest_gc_seq['gc']

        # ==========================
        # GC DISTRIBUTION
        # ==========================

        low_gc = 0
        medium_gc = 0
        high_gc = 0

        for item in sorted_data:

            if item['gc'] < 40:

                low_gc += 1

            elif item['gc'] <= 60:

                medium_gc += 1

            else:

                high_gc += 1

        # ==========================
        # NUCLEOTIDE PERCENTAGE
        # ==========================

        total_bases = sum(
            total_freq.values()
        )

        nucleotide_percent = {
            'A': 0,
            'C': 0,
            'G': 0,
            'T': 0
        }

        if total_bases > 0:

            nucleotide_percent = {

                nuc:
                (count / total_bases) * 100

                for nuc, count
                in total_freq.items()

            }

        # ==========================
        # GRAPH
        # ==========================

        plot_url = create_plot(
            sorted_data
        )

        # ==========================
        # CSV EXPORT
        # ==========================

        csv_path = os.path.join(
            'results',
            'analysis_results.csv'
        )

        with open(
            csv_path,
            'w',
            newline='',
            encoding='utf-8'
        ) as f:

            writer = csv.writer(f)

            writer.writerow([

                'File',
                'Sequence ID',
                'Length',
                'GC Content (%)',
                'A',
                'C',
                'G',
                'T'

            ])

            for item in sorted_data:

                writer.writerow([

                    item['file'],
                    item['id'],
                    item['length'],
                    f"{item['gc']:.2f}",
                    item['freq']['A'],
                    item['freq']['C'],
                    item['freq']['G'],
                    item['freq']['T']

                ])

        return render_template(

            'index.html',

            top_3=top_3,
            top_10=top_10,

            total_sequences=total_sequences,

            avg_gc=avg_gc,
            avg_length=avg_length,

            max_gc=max_gc,
            min_gc=min_gc,

            highest_gc_seq=highest_gc_seq,
            lowest_gc_seq=lowest_gc_seq,

            low_gc=low_gc,
            medium_gc=medium_gc,
            high_gc=high_gc,

            total_freq=total_freq,
            nucleotide_percent=nucleotide_percent,

            all_data=sorted_data,

            plot_url=plot_url,
            csv_ready=True

        )

    return render_template(
        'index.html',
        top_3=None
    )


@app.route('/download')
def download():

    return send_file(
        'results/analysis_results.csv',
        as_attachment=True
    )


if __name__ == '__main__':
    app.run(debug=True)
