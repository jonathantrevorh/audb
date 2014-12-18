.echo on

create virtual table songs using audb_tree;
insert into songs values('1', '"./Music/Aoki - Boneless.mp3"', 'Boneless', 'Steve Aoki');
select * from songs where song = '"./Music/Aoki2.wav"';
