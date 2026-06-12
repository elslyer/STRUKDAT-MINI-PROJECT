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

    top_plot = data[:20]

    ids = [item['id'][:10] for item in top_plot]
    gc_values = [item['gc'] for item in top_plot]

    plt.figure(figsize=(12, 6))
    plt.bar(ids, gc_values)

    plt.xlabel("Sequence ID")
    plt.ylabel("GC Content (%)")
    plt.title("GC Content per Sequence")
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

        top_3 = sorted_data[:3]

        total_sequences = len(sorted_data)

        if total_sequences > 0:

            avg_gc = sum(
                item['gc']
                for item in sorted_data
            ) / total_sequences

            max_gc = max(
                item['gc']
                for item in sorted_data
            )

            min_gc = min(
                item['gc']
                for item in sorted_data
            )

        else:

            avg_gc = 0
            max_gc = 0
            min_gc = 0

        plot_url = create_plot(sorted_data)

        csv_path = os.path.join(
            'results',
            'analysis_results.csv'
        )

        with open(
            csv_path,
            'w',
            newline=''
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
            plot_url=plot_url,
            csv_ready=True,
            total_sequences=total_sequences,
            avg_gc=avg_gc,
            max_gc=max_gc,
            min_gc=min_gc,
            total_freq=total_freq
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