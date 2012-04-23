        #!/usr/bin/perl -w
        use strict;
        use IPC::Shareable;
        my $glue = 'data';
        my %options = (
            create    => 'yes',
            exclusive => 0,
            mode      => 0666,
            destroy   => 'yes',
        );
        my %colours;
        tie %colours, 'IPC::Shareable', $glue, { %options } or
            die "server: tie failed\n";
        %colours = (
            red => [
                'fire truck',
                'leaves in the fall',
            ],
            blue => [
                'sky',
                'police cars',
            ],
        );


unless (my $child = fork)
{
        my %options = (
            create    => 0,
            exclusive => 0,
            mode      => 0644,
            destroy   => 0,
        );
        my %colours;
        tie %colours, 'IPC::Shareable', $glue, { %options } or die "tie failed\n";
        $colours{black} = ['night', 'ink'];
        exit;	
}

        print "server: here are all my colours:\n";
        foreach my $c (keys %colours) {
            print "server: these are $c: ",
                join(', ', @{$colours{$c}}), "\n";
        }
        exit;
