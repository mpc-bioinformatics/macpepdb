import re
import psycopg2

from ..models.peptide import Peptide

class Statistics:
    # Can be overriden in subclasses to change the peptide tables
    peptide_class=Peptide

    @classmethod
    def estimate_peptide_partition_utilizations(cls, database_cursor) -> list:
        """
        Estimates the parititon utilization. Counting peptides need to much time, so we can estimate it with help of the pg_class view.
        The estimation changes a bit after each VACUUM or ANALYZE but is very fast.
        @param database_cursor
        @return List of tupels [(parition_count, partition_name), ...]
        """
        database_cursor.execute(f"SELECT relname, reltuples::BIGINT FROM pg_class WHERE relname SIMILAR TO '{cls.peptide_class.TABLE_NAME}_[0-9]{{3}}';")
        return database_cursor.fetchall()

    @classmethod
    def estimate_peptide_count(cls, database_cursor) -> int:
        """
        Estimates the peptide count, by using the partition utilizaton.
        The estimation changes a bit after each VACUUM or ANALYZE but is very fast.
        @param database_cursor
        @return Estimated peptide count
        """
        sum = 0
        for partition_estimation in cls.estimate_peptide_partition_utilizations(database_cursor):
            sum += partition_estimation[1]
        return sum

    @classmethod
    def get_partition_boundaries(cls, database_cursor):
        """
        This return the peptide partition boundaries.
        @param database_cursor
        @return List of tupel [(partition_name, from, to), ...]
        """
        database_cursor.execute(f"select pg_class.relname, pg_get_expr(pg_class.relpartbound, pg_class.oid, true) from pg_class where relname SIMILAR TO '{cls.peptide_class.TABLE_NAME}_[0-9]{{3}}';")
        rows = database_cursor.fetchall()
        num_regex = re.compile(r"\d+")
        partition_boundaries = []
        for row in rows:
            matches = re.findall(num_regex, row[1])
            partition_boundaries.append((row[0], int(matches[0]), int(matches[1])))
        partition_boundaries.sort(key=lambda partition_boundary: partition_boundary[1])
        return partition_boundaries

    @classmethod
    def get_partition_boundaries_from_command_line(cls, args):
        database_connection = psycopg2.connect(args.database_url)
        with database_connection:
            with database_connection.cursor() as database_cursor:
                partition_boundaries = cls.get_partition_boundaries(database_cursor)
                if not args.csv:
                    for partition in partition_boundaries:
                        print(f"{partition[0]:<15}\t{partition[1]:<15}\t{partition[2]:<15}")
                else:
                    print("\"parition name\", \"lower boundary\", \"upper boundary\"")
                    for partition in partition_boundaries:
                        print(f"\"{partition[0]}\", {partition[1]}, {partition[2]}")
        database_connection.close()



    @classmethod
    def __comand_line_arguments_for_partition_boundaries(cls, subparsers):
        parser = subparsers.add_parser('peptide-partition-boundaries', help="Prints the used peptide partition boundaries.")
        parser.add_argument("--csv", action="store_const", const="True", help="Output is in csv format.", default=False)
        parser.set_defaults(func=cls.get_partition_boundaries_from_command_line)

    @classmethod
    def estimate_partition_usage_from_comand_line(cls, args):
        database_connection = psycopg2.connect(args.database_url)
        with database_connection:
            with database_connection.cursor() as database_cursor:
                parition_estimations = cls.estimate_peptide_partition_utilizations(database_cursor)
                if not args.csv:
                    for parition_estimation in parition_estimations:
                        print(f"{parition_estimation[0]:<15}\t{parition_estimation[1]:<15}")
                else:
                    print("\"parition name\", \"count\"")
                    for parition_estimation in parition_estimations:
                        print(f"\"{parition_estimation[0]}\", {parition_estimation[1]}")
        database_connection.close()

    @classmethod
    def __comand_line_arguments_for_patition_usage(cls, subparsers):
        parser = subparsers.add_parser('estimate-peptide-partition-utilizations', help="Prints estimation for peptide parition utilization.")
        parser.add_argument("--csv", action="store_const", const="True", help="Output is in csv format.", default=False)
        parser.set_defaults(func=cls.estimate_partition_usage_from_comand_line)

    @classmethod
    def count_peptides_from_command_line(cls, args):
        database_connection = psycopg2.connect(args.database_url)
        with database_connection:
            with database_connection.cursor() as database_cursor:
                print(cls.estimate_peptide_count(database_cursor))
        database_connection.close()

    @classmethod
    def __comand_line_arguments_peptide_counts(cls, subparsers):
        parser = subparsers.add_parser('estimate-peptide-count', help="Estimates the peptide count")
        parser.set_defaults(func=cls.count_peptides_from_command_line)

    @classmethod
    def command_line_arguments(cls, subparsers):
        parser = subparsers.add_parser('statistics', help="Gather some database statistics")
        parser.add_argument("--database-url", "-d", type=str, required=True, help="Database URL for postgres, e.g. postgres://user:password@server/database")
        subsubparsers = parser.add_subparsers()
        cls.__comand_line_arguments_for_patition_usage(subsubparsers)
        cls.__comand_line_arguments_peptide_counts(subsubparsers)
        cls.__comand_line_arguments_for_partition_boundaries(subsubparsers)