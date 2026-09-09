import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0005_add_ml_models'),
    ]

    operations = [
        migrations.CreateModel(
            name='GEOExpressionRecord',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('gene_symbol', models.CharField(db_index=True, max_length=50)),
                ('value', models.FloatField()),
                ('zscore', models.FloatField(blank=True, null=True)),
                ('geo_sample', models.ForeignKey(
                    on_delete=django.db.models.deletion.CASCADE,
                    related_name='expression_records',
                    to='core.geosample',
                )),
                ('hgnc_gene', models.ForeignKey(
                    blank=True,
                    null=True,
                    on_delete=django.db.models.deletion.CASCADE,
                    related_name='geo_expression_records',
                    to='core.hgncgene',
                )),
            ],
            options={
                'db_table': 'geo_expression_record',
            },
        ),
        migrations.AddIndex(
            model_name='geoexpressionrecord',
            index=models.Index(fields=['gene_symbol'], name='geo_expr_symbol_idx'),
        ),
        migrations.AddIndex(
            model_name='geoexpressionrecord',
            index=models.Index(fields=['geo_sample', 'gene_symbol'], name='geo_expr_sample_gene_idx'),
        ),
        migrations.AddConstraint(
            model_name='geoexpressionrecord',
            constraint=models.UniqueConstraint(
                fields=['geo_sample', 'gene_symbol'],
                name='uq_geo_expr_sample_gene',
            ),
        ),
    ]
